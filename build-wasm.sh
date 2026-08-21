#!/bin/bash
# Build Corset (BAM-free) as a DCP-ready WASM module.
# Prereq: emsdk activated; run from the Corset repo with corset-nobam.patch applied.
set -e

em++ -O2 -std=c++11 -DUNORDEREDMAP -DVERSION=1.09 -I. \
  corset-nobam.cc Cluster.cc MakeClusters.cc Transcript.cc Read.cc \
  -s MODULARIZE=1 -s EXPORT_NAME=CorsetModule \
  -s SINGLE_FILE=1 \
  -s ALLOW_MEMORY_GROWTH=1 \
  -s INVOKE_RUN=0 \
  -s EXPORTED_RUNTIME_METHODS='["callMain","FS"]' \
  -s EXIT_RUNTIME=1 \
  -o corset-module.js

# Patch node:-scheme requires so DCP's dependency bundler can resolve them
# (per https://docs.dcp.dev/advanced/wasm-modules.html)
node -e '
const fs = require("fs");
let d = fs.readFileSync("corset-module.js", "utf8");
d = d.replaceAll("require(\"node:fs\")", "require(\"fs\")")
     .replaceAll("require(\"node:crypto\")", "require(\"crypto\")");
fs.writeFileSync("corset-module.js", d);
'
echo "Built corset-module.js"
