// corset-wasm.js -- friendly wrapper: eq_classes text in, clusters/counts out.
// Each call gets a fresh module instance (Corset's main() assumes a clean run).
const CorsetModule = require('./corset-module.js');

/**
 * Run Corset on salmon-style eq_classes content.
 * @param {string[]} eqClassFiles - one string of eq_classes.txt content per sample
 * @param {string[]} extraArgs    - e.g. ['-d','0.3'] for a distance threshold
 * @returns {{clusters: string, counts: string, log: string}}
 */
async function runCorset(eqClassFiles, extraArgs = []) {
  let log = '';
  const Module = await CorsetModule({
    print: (s) => (log += s + '\n'),
    printErr: (s) => (log += '[err] ' + s + '\n'),
  });
  const names = eqClassFiles.map((_, i) => `/sample${i}.eq`);
  names.forEach((n, i) => Module.FS.writeFile(n, eqClassFiles[i]));
  Module.callMain(['-i', 'salmon_eq_classes', ...extraArgs, ...names]);
  return {
    clusters: Module.FS.readFile('clusters.txt', { encoding: 'utf8' }),
    counts: Module.FS.readFile('counts.txt', { encoding: 'utf8' }),
    log,
  };
}

module.exports = { runCorset };
