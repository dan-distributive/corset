/**
 * @file        bootstrap-corset-job.js
 * @description Parameter-sweep + bootstrap cluster-stability harness for the
 *              WASM port of Corset (see corset-nobam.patch, build-wasm.sh,
 *              corset-wasm.js). Answers a real question a Corset user has to
 *              answer by hand today: "which -d distance threshold gives
 *              clusters that are actually reproducible, not just an
 *              artifact of this one read sample?"
 *
 *              For each candidate threshold, dispatches many independent
 *              Poisson-bootstrap resamples of the input equivalence classes
 *              (classic non-parametric bootstrap over read counts -- same
 *              idea as edgeR/DESeq2 count bootstrapping) plus one
 *              unresampled "point estimate" run, as a single DCP fleet. Each
 *              worker resamples its own replicate and runs actual Corset
 *              (via corset-wasm.js, published as the 'corset-test' package --
 *              see package/package.dcp) -- the resampling has to happen on
 *              the worker, not the client, or you're back to shipping
 *              thousands of resampled eq_classes files over the wire.
 *              The client then scores each threshold by mean Adjusted Rand
 *              Index (ARI) between each replicate's clustering and that
 *              threshold's point estimate: near 1.0 means bootstrap noise
 *              barely moves the clusters (stable); near 0 means the
 *              clustering is not trustworthy at that threshold.
 *
 *              Modeled on job-templates/node/job-template.js's structure
 *              (wallet/identity/compute.for/job.public/events) and
 *              ffmpeg-dcp/app.js's dispatchJob (work function does the
 *              worker-side heavy lifting; client aggregates results as they
 *              stream in). Prerequisite: publish package/corset-wasm.js as
 *              'corset-test' via dcp-util's `publish package
 *              /absolute/path/to/package/package.dcp` -- see
 *              [[dcp_wasm_packaging_howto]] for the publish gotchas
 *              (job.requires() auto-discovery only goes ~2 deep, but this
 *              bundle is a single flattened file so that doesn't bite here).
 *
 * @usage       node bootstrap-corset-job.js --apiKey=0x<revocable key> [sample1.eq sample2.eq ...] \
 *                  [--thresholds=0.2,0.3,0.4] [--bootstraps=20] [--seed=1] [--computeGroup=key,secret]
 *              With no sample files given, runs against small synthetic
 *              eq_classes data so the harness itself can be smoke-tested
 *              without a real salmon run.
 *
 *              eq_classes file format (salmon's, same as corset -i
 *              salmon_eq_classes expects):
 *                <ntranscripts>
 *                <neqclasses>
 *                <transcript id> * ntranscripts
 *                <eq_size> <trans_pos>*eq_size <weight>   * neqclasses
 */
'use strict';

// ---------------------------------------------------------------------------
// eq_classes parsing/serialization + seeded Poisson bootstrap. Duplicated
// (not shared via require) between here and workFunction below on purpose --
// DCP work functions are re-evaluated in an isolated worker context from
// their own toString(), so helpers used inside workFunction must be nested
// inside it (or pulled in via job.requires()), not imported from module
// scope. See [[dcp_wasm_packaging_howto]].
// ---------------------------------------------------------------------------

function parseEqClasses(text) {
  const lines = text.split('\n');
  let i = 0;
  const ntranscripts = parseInt(lines[i++], 10);
  const neqclasses = parseInt(lines[i++], 10);
  const transIds = [];
  for (let t = 0; t < ntranscripts; t++) transIds.push(lines[i++]);
  const classes = [];
  for (let e = 0; e < neqclasses; e++) {
    const parts = lines[i++].trim().split(/\s+/).map(Number);
    const eqSize = parts[0];
    const transPos = parts.slice(1, 1 + eqSize);
    const weight = parts[1 + eqSize];
    classes.push({ transPos, weight });
  }
  return { transIds, classes };
}

function serializeEqClasses({ transIds, classes }) {
  const lines = [String(transIds.length), String(classes.length), ...transIds];
  for (const c of classes) {
    lines.push([c.transPos.length, ...c.transPos, c.weight].join(' '));
  }
  return lines.join('\n') + '\n';
}

// mulberry32: tiny seeded PRNG, deterministic per (thresholdIndex, seed,
// sampleIndex) so a given replicate is exactly reproducible on re-dispatch.
function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// Knuth's algorithm -- fine for the read-count magnitudes eq_classes weights
// take (not appropriate for lambda in the many-thousands, where this loop
// gets slow; swap for a normal-approximation sampler if that matters for a
// given dataset).
function poissonSample(lambda, rng) {
  if (lambda <= 0) return 0;
  const L = Math.exp(-lambda);
  let k = 0, p = 1;
  do { k++; p *= rng(); } while (p > L);
  return k - 1;
}

// Poisson bootstrap: each equivalence class's read count is independently
// resampled as Poisson(original count). Standard approximation to the
// multinomial read-resampling bootstrap for count data at this scale, and
// embarrassingly parallel across classes (no shared resampling state
// needed, unlike a true multinomial draw).
function poissonResample(parsed, rng) {
  return {
    transIds: parsed.transIds,
    classes: parsed.classes.map((c) => ({ transPos: c.transPos, weight: poissonSample(c.weight, rng) })),
  };
}

function parseClustersText(clustersTxt) {
  const map = new Map();
  for (const line of clustersTxt.split('\n')) {
    if (!line.trim()) continue;
    const [transcript, cluster] = line.split('\t');
    map.set(transcript, cluster);
  }
  return map;
}

// Adjusted Rand Index between two clusterings of the same item set (here:
// transcripts). 1.0 = identical partitions, ~0 = agreement no better than
// chance. Standard bootstrap cluster-stability metric.
function adjustedRandIndex(mapA, mapB) {
  const items = [...mapA.keys()].filter((k) => mapB.has(k));
  const n = items.length;
  if (n < 2) return 1;

  const countsA = new Map(); // clusterLabel -> count
  const countsB = new Map();
  const contingency = new Map(); // "labelA|labelB" -> count
  for (const item of items) {
    const a = mapA.get(item), b = mapB.get(item);
    countsA.set(a, (countsA.get(a) || 0) + 1);
    countsB.set(b, (countsB.get(b) || 0) + 1);
    const key = a + '|' + b;
    contingency.set(key, (contingency.get(key) || 0) + 1);
  }
  const choose2 = (x) => (x * (x - 1)) / 2;

  let sumComb = 0;
  for (const c of contingency.values()) sumComb += choose2(c);
  let sumA = 0;
  for (const c of countsA.values()) sumA += choose2(c);
  let sumB = 0;
  for (const c of countsB.values()) sumB += choose2(c);
  const totalComb = choose2(n);

  const expectedIndex = (sumA * sumB) / totalComb;
  const maxIndex = 0.5 * (sumA + sumB);
  if (maxIndex === expectedIndex) return 1; // degenerate (e.g. every item its own cluster)
  return (sumComb - expectedIndex) / (maxIndex - expectedIndex);
}

// ---------------------------------------------------------------------------
// Synthetic demo data, used only when no real eq_classes files are passed on
// the command line -- lets the harness itself be smoke-tested end-to-end.
// ---------------------------------------------------------------------------

function makeSyntheticSample(baseWeight, jitter, rngSeed) {
  const rng = mulberry32(rngSeed);
  // 6 transcripts across 2 well-separated true clusters {t0,t1,t2} / {t3,t4,t5},
  // plus cross-cluster noise so threshold choice actually matters.
  const transIds = ['t0', 't1', 't2', 't3', 't4', 't5'];
  const classes = [
    { transPos: [0, 1], weight: Math.round(baseWeight * (1 + jitter * (rng() - 0.5))) },
    { transPos: [1, 2], weight: Math.round(baseWeight * (1 + jitter * (rng() - 0.5))) },
    { transPos: [0, 2], weight: Math.round(baseWeight * 0.6 * (1 + jitter * (rng() - 0.5))) },
    { transPos: [3, 4], weight: Math.round(baseWeight * (1 + jitter * (rng() - 0.5))) },
    { transPos: [4, 5], weight: Math.round(baseWeight * (1 + jitter * (rng() - 0.5))) },
    { transPos: [3, 5], weight: Math.round(baseWeight * 0.6 * (1 + jitter * (rng() - 0.5))) },
    { transPos: [2, 3], weight: Math.round(baseWeight * 0.15 * (1 + jitter * (rng() - 0.5))) }, // cross-cluster noise
  ];
  return serializeEqClasses({ transIds, classes });
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

async function main() {
  const identity = require('dcp/identity');
  const compute = require('dcp/compute');
  const fs = require('fs');

  function getFlag(name, def) {
    const pfx = `--${name}=`;
    const hit = process.argv.find((a) => a.startsWith(pfx));
    return hit ? hit.slice(pfx.length) : def;
  }

  const thresholds = getFlag('thresholds', '0.2,0.3,0.4').split(',').map(Number);
  const nBootstrap = parseInt(getFlag('bootstraps', '20'), 10);
  const baseSeed = parseInt(getFlag('seed', '1'), 10);
  const apiKey = getFlag('apiKey');
  const computeGroupFlag = getFlag('computeGroup');
  const sampleFiles = process.argv.slice(2).filter((a) => !a.startsWith('--'));

  if (!apiKey) { console.error('ERROR: --apiKey=0x... is required to launch a real DCP job'); process.exit(1); }

  /* INPUT DATA */
  let eqClassFiles;
  if (sampleFiles.length > 0) {
    eqClassFiles = sampleFiles.map((f) => fs.readFileSync(f, 'utf8'));
  } else {
    console.log('No sample files given -- using synthetic demo data (2 samples, 6 transcripts, 2 true clusters).');
    eqClassFiles = [makeSyntheticSample(40, 0.3, 1001), makeSyntheticSample(35, 0.3, 2002)];
  }
  const eqClassFilesJson = JSON.stringify(eqClassFiles);

  /* DCP API KEY -- a raw private key used directly as both identity and
     payer, the same pattern proven in bootstrap-realization-launcher.js's
     real-dispatch branch (identity.set(apiKey) + job.exec(slicePrice), no
     separate wallet.get('default')/setPaymentAccountKeystore needed). */
  await identity.set(apiKey);

  /* INPUT SET: cartesian product of thresholds x (point-estimate + bootstrap replicates) */
  const inputSet = [];
  for (const threshold of thresholds) {
    inputSet.push({ threshold, seed: null, isPointEstimate: true });
    for (let b = 0; b < nBootstrap; b++) {
      inputSet.push({ threshold, seed: baseSeed * 1000 + b, isPointEstimate: false });
    }
  }

  /* WORK FUNCTION -- runs on the DCP worker */
  async function workFunction(unit, eqClassFilesJsonArg) {
    progress();

    function parseEqClasses(text) {
      const lines = text.split('\n');
      let i = 0;
      const ntranscripts = parseInt(lines[i++], 10);
      const neqclasses = parseInt(lines[i++], 10);
      const transIds = [];
      for (let t = 0; t < ntranscripts; t++) transIds.push(lines[i++]);
      const classes = [];
      for (let e = 0; e < neqclasses; e++) {
        const parts = lines[i++].trim().split(/\s+/).map(Number);
        const eqSize = parts[0];
        const transPos = parts.slice(1, 1 + eqSize);
        const weight = parts[1 + eqSize];
        classes.push({ transPos, weight });
      }
      return { transIds, classes };
    }
    function serializeEqClasses({ transIds, classes }) {
      const lines = [String(transIds.length), String(classes.length), ...transIds];
      for (const c of classes) lines.push([c.transPos.length, ...c.transPos, c.weight].join(' '));
      return lines.join('\n') + '\n';
    }
    function mulberry32(seed) {
      let a = seed >>> 0;
      return function () {
        a |= 0; a = (a + 0x6d2b79f5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
      };
    }
    function poissonSample(lambda, rng) {
      if (lambda <= 0) return 0;
      const L = Math.exp(-lambda);
      let k = 0, p = 1;
      do { k++; p *= rng(); } while (p > L);
      return k - 1;
    }
    function poissonResample(parsed, rng) {
      return {
        transIds: parsed.transIds,
        classes: parsed.classes.map((c) => ({ transPos: c.transPos, weight: poissonSample(c.weight, rng) })),
      };
    }

    // Resolved via job.requires(['corset-test/corset-wasm.js']) below.
    const { runCorset } = require('corset-wasm.js');

    const eqClassFilesArg = JSON.parse(eqClassFilesJsonArg);
    const { threshold, seed } = unit;

    let samples = eqClassFilesArg;
    if (seed !== null) {
      // Independent resample per sample, seeded off (replicate seed, sample
      // index) so samples in a replicate aren't spuriously correlated.
      samples = eqClassFilesArg.map((text, sampleIndex) => {
        const rng = mulberry32(seed * 100003 + sampleIndex);
        return serializeEqClasses(poissonResample(parseEqClasses(text), rng));
      });
    }

    const { clusters } = await runCorset(samples, ['-d', String(threshold)]);
    return { threshold, seed, clusters };
  }

  /* JOB */
  const job = compute.for(inputSet, workFunction, [eqClassFilesJson]);
  job.requires(['corset-test/corset-wasm.js']);

  /* COMPUTE GROUPS */
  job.computeGroups = (computeGroupFlag || 'public').split(':').map((g) => g.trim()).filter(Boolean).map((group) => {
    const [joinKey, joinSecret] = group.split(',').map((s) => s.trim());
    return joinSecret ? { joinKey, joinSecret } : { joinKey };
  });

  /* PUBLIC INFO */
  job.public = {
    name: '🧬 Corset bootstrap cluster-stability sweep',
    description: `${thresholds.length} thresholds x (1 point estimate + ${nBootstrap} bootstrap replicates)`,
  };

  /* EVENTS: aggregate ARI per threshold as results stream in */
  const pointEstimates = new Map();  // threshold -> Map<transcript, cluster>
  const replicateResults = new Map(); // threshold -> [Map<transcript, cluster>, ...]
  for (const t of thresholds) replicateResults.set(t, []);

  function handleResult({ threshold, seed, clusters }) {
    const clusterMap = parseClustersText(clusters);
    if (seed === null) {
      pointEstimates.set(threshold, clusterMap);
      console.log(`  [threshold ${threshold}] point estimate received (${clusterMap.size} transcripts)`);
    } else {
      replicateResults.get(threshold).push(clusterMap);
      console.log(`  [threshold ${threshold}] replicate ${replicateResults.get(threshold).length}/${nBootstrap} received`);
    }
  }

  job.on('readystatechange', (ev) => console.log(`Ready state: ${ev}`));
  job.on('accepted', () => console.log(`  Job id: ${job.id}\n  Awaiting results...`));
  job.on('result', (ev) => handleResult(ev.result));
  job.on('error', (error) => console.error('  Job error:', error));
  job.on('nofunds', (ev) => console.log(ev));

  /* EXECUTION -- apiKey is both identity and payer, so no separate
     wallet/payment-keystore step; job.exec(slicePrice) pays from it directly. */
  const slicePrice = compute.marketValue(1.0);
  await job.exec(slicePrice);

  /* REPORT: mean ARI (bootstrap replicate vs. point estimate) per threshold */
  console.log('\n=== Cluster stability by threshold (mean ARI, bootstrap vs. point estimate) ===\n');
  for (const threshold of thresholds) {
    const point = pointEstimates.get(threshold);
    const reps = replicateResults.get(threshold);
    if (!point || reps.length === 0) {
      console.log(`  d=${threshold}: no data`);
      continue;
    }
    const aris = reps.map((rep) => adjustedRandIndex(point, rep));
    const mean = aris.reduce((a, b) => a + b, 0) / aris.length;
    const sorted = [...aris].sort((a, b) => a - b);
    console.log(`  d=${threshold}: mean ARI=${mean.toFixed(3)}  median=${sorted[Math.floor(sorted.length / 2)].toFixed(3)}  min=${sorted[0].toFixed(3)}  (n=${aris.length})`);
  }
  console.log('\nHigher mean ARI = clustering at that threshold is more robust to resampling noise in the input reads.');
}

require('dcp-client')
  .init('https://scheduler.distributed.computer')
  .then(main)
  .catch(console.error);
