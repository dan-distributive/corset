# corset (WASM / DCP port)

A [WebAssembly](https://webassembly.org/) build of [Oshlack/Corset](https://github.com/Oshlack/Corset) —
the de novo transcriptome clustering/counting tool — packaged to run as a
[Distributed Compute Protocol (DCP)](https://dcp.dev) job, plus a bootstrap
cluster-stability harness built on top of it.

Corset takes multi-mapped read alignments (or salmon equivalence classes)
from a de novo assembled transcriptome and clusters contigs that are likely
to be the same gene, then produces per-cluster read counts for downstream
differential expression analysis. This port only supports the
`salmon_eq_classes` input path — see [Why BAM support was dropped](#why-bam-support-was-dropped)
below.

## How this was built

This module was produced with [Claude Code](https://claude.com/claude-code),
following the WASM-packaging recipe in
[docs.dcp.dev/advanced/wasm-modules.html](https://docs.dcp.dev/advanced/wasm-modules.html).
Every step below — including finding and fixing two real bugs along the
way — was done and verified end-to-end (compiled, smoke-tested locally,
then dispatched as a real job against the public DCP network) in the same
session, not just written from the docs and left unverified.

1. **Strip the BAM/samtools dependency.** Corset's BAM-reading path
   (`read_bam_file()` in `corset.cc`) requires `libbam`/`libsam`, C libraries
   that don't cross-compile to WASM cleanly and aren't needed for the
   `salmon_eq_classes` input path this port targets. `corset-nobam.patch`
   guards the BAM code out behind `#ifdef HAVE_BAM` (left undefined at build
   time) and adds a stub that exits with a clear error if `-i bam` is ever
   selected. It also adds a missing `#include <cstdint>` to `Read.h`, needed
   under Emscripten's newer libc++ headers.

2. **Compile to WASM** (`build-wasm.sh`), following docs.dcp.dev's recipe:
   `-s MODULARIZE=1` + `-s EXPORT_NAME=CorsetModule` so the build is a
   self-contained factory function, `-s SINGLE_FILE=1` so the `.wasm` binary
   is base64-embedded directly in the `.js` (no second file to ship/fetch —
   DCP workers can't `fetch()` arbitrary URLs), and `-s INVOKE_RUN=0` so the
   module doesn't auto-run `main()` on instantiation; the wrapper drives
   Corset's `main()` explicitly via `callMain()` after staging input files
   in Emscripten's in-memory filesystem (MEMFS).

   **Bug found here:** the first draft of this script used `emcc` (the C
   compiler driver) instead of `em++` to compile/link Corset's `.cc` files.
   `emcc` doesn't link libc++, so the link step failed with dozens of
   undefined `std::__2::*` symbols (`ios_base`, `ostream`, `locale`, `operator
   new`/`delete`, ...). One-word fix: `em++`.

3. **Wrap it for a normal JS call shape** (`corset-wasm.js`): `runCorset(eqClassFiles, extraArgs)`
   writes each sample's eq_classes content into MEMFS, calls `callMain(['-i',
   'salmon_eq_classes', ...extraArgs, ...files])`, and returns
   `{ clusters, counts, log }` read back out of MEMFS — no Emscripten,
   `ccall`/MEMFS, or Corset CLI knowledge required of a caller.

4. **Bundle it for DCP's package manager** (`corset-package/`):
   DCP's `job.requires()` can pull in a published, pre-built package so a
   work function can just `require('corset-wasm.js')` on any worker without
   shipping the ~300 KB module as job arguments on every dispatch.
   `build-bravojs-bundle.js` flattens `corset-module.js` +
   `corset-wasm.js` into one file wrapped in
   [bravojs](https://github.com/Distributive-Network/bravojs)'s
   `module.declare(...)` format, which is what the package publisher expects
   (local/unpublished `job.requires()` tolerates plain CommonJS because it
   goes through a real webpack build at deploy time; publishing does not).

   **Bug found here:** the manifest's version was written as `1.09.0`,
   which looks like "Corset 1.09" but isn't valid semver — leading zeros
   aren't allowed in numeric identifiers, so `publish` rejected it with
   `EBADPARAM`. Fixed to `1.9.0`.

   Published via `dcp-util`'s publish tool:
   ```
   node node_modules/dcp-util/bin/publish package /absolute/path/to/corset-package/package.dcp
   ```
   (path must be absolute — it resolves against the shell's cwd, not the
   manifest's own location. No `--apiKey` needed for publishing itself,
   only for dispatching jobs.)

### Why BAM support was dropped

DCP workers run in a sandboxed environment and receive job input as data,
not as files on a real filesystem `libbam` could open — and multi-mapped
BAM files for a real dataset are typically much larger than the
`salmon_eq_classes` files that already encode the same equivalence-class
information Corset actually clusters on. Rather than trying to cross-compile
`libbam`/`zlib` and design a way to stream large BAMs into MEMFS, this port
targets the `-i salmon_eq_classes` path directly, which is already the
recommended input format for de novo assemblies without a reference genome.

### License

Corset itself is GPL-3.0 ([COPYING](https://github.com/Oshlack/Corset/blob/master/COPYING)
in the upstream repo). This is an unmodified-algorithm, patched-for-portability
build (see `corset-nobam.patch`) distributed as source, which satisfies
GPL's requirements; it isn't folded into any proprietary code.

## `bootstrap-corset-job.js`: a real use for it

Corset's `-d` distance threshold controls how aggressively contigs get
merged into the same cluster — but there's no built-in way to tell whether
a given threshold's clusters are actually reproducible, or just an artifact
of one particular read sample. This harness answers that:

For each candidate threshold, it dispatches a single DCP fleet containing
one unresampled "point estimate" run plus many independent
[Poisson-bootstrap](https://en.wikipedia.org/wiki/Bootstrapping_(statistics))
resamples of the input equivalence classes (the same idea edgeR/DESeq2 use
for count-based bootstrapping). Each worker resamples its own replicate
*and* runs actual Corset via the published `corset-wasm.js` package — the
resampling has to happen on the worker, not the client, or you're back to
shipping thousands of resampled files over the wire. The client then scores
each threshold by mean [Adjusted Rand Index](https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index)
(ARI) between each replicate's clustering and that threshold's point
estimate — near 1.0 means bootstrap noise barely moves the clusters
(stable and trustworthy); lower means that threshold's clusters shouldn't
be trusted at face value.

With no sample files given, it runs against small synthetic eq_classes data
(2 samples, 6 transcripts, 2 true clusters) so the harness itself can be
smoke-tested without a real salmon run.

### Usage

```
node bootstrap-corset-job.js [sample1.eq sample2.eq ...] \
  --apiKey=0x<your DCP identity key> \
  [--thresholds=0.2,0.3,0.4] [--bootstraps=20] [--seed=1] [--computeGroup=key,secret]
```

Example (revocable test key and private compute group obfuscated):

```
dandesjardins@Dans-MacBook-Air-2 corset % node bootstrap-corset-job.js \
  --thresholds=0.2,0.3,0.4 \
  --bootstraps=20 \
  --seed=1 \
  --apiKey=0x<redacted> \
  --computeGroup=<redacted>,<redacted>
```

Output is a per-threshold stability report:

```
=== Cluster stability by threshold (mean ARI, bootstrap vs. point estimate) ===

  d=0.2: mean ARI=1.000  median=1.000  min=1.000  (n=20)
  d=0.3: mean ARI=0.933  median=0.950  min=0.812  (n=20)
  d=0.4: mean ARI=0.788  median=0.801  min=0.640  (n=20)

Higher mean ARI = clustering at that threshold is more robust to resampling noise in the input reads.
```
