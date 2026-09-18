#!/bin/bash
# matchall.sh <out-dir> <sample>... : deposited reference vs solved res, one summary line each
export PYTHONHOME='C:\Users\florian\AppData\Local\Python\pythoncore-3.12-64' LIBTBX_BUILD='D:\devel\rundir-py3\cctbx\cctbx_build' PYTHONPATH='D:\devel\rundir-py3\cctbx\cctbx_sources;D:\devel\rundir-py3\cctbx\cctbx_sources\boost_adaptbx;D:\devel\rundir-py3\cctbx\cctbx_build\lib;D:\devel\rundir-py3\util\pyUtil\Lib\site-packages'
PY=/c/Users/florian/AppData/Local/Python/pythoncore-3.12-64/python.exe
out=$1; shift
for s in "$@"; do
  ref=""
  for d in /d/devel/olex2-samples/$s /d/devel/rundir-test/sample_data/$s; do
    for e in ins res cif; do f=$(ls $d/*.$e 2>/dev/null | head -1); [ -n "$f" ] && { ref=$f; break 2; }; done
  done
  res=$(ls $out/scratch/$s/*.res 2>/dev/null | head -1)
  echo "## $s ref=$(basename "$ref") res=$(basename "$res")"
  [ -n "$ref" ] && [ -n "$res" ] && "$PY" "$(dirname "$0")/match.py" "$ref" "$res" 2>&1 | grep -E "shift|<--|Error|error" | head -20
done
