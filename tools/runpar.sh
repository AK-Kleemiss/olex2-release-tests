#!/bin/bash
# runpar.sh <pass-name> <ids-file> <nworkers> : Auto-Solve release cases spread over rundir-w1..wN in parallel, then match + summary
S="C:/Users/florian/AppData/Local/Temp/claude/E--Dropbox-Obsidian-LLM-Infos/ab19d726-5706-4c7b-994a-827815740ab3/scratchpad"
name=$1; ids=$2; nw=$3
split -n l/$nw -d "$ids" "$S/.chunk_${name}_"
k=1
for c in "$S"/.chunk_${name}_*; do
  rd="${RD:-D:/devel/rundir-w}$k"; out="D:/devel/olex2-release-out/${name}_w$k"
  rm -f "${rd/D:/\/d}/runonce.release_tests.txm"
  list=$(tr '\n' ',' < "$c" | sed 's/,$//')
  pwsh -NoProfile -Command "\$s='$list'; \$env:OLEX2_TEST_AUTOSOLVE_SAMPLES=\$s; \$cases=(\$s.Split(',') | ForEach-Object { \"autosolve_\$(\$_.ToLower())\" }) -join ','; Set-Location D:\devel\olex2-release-tests; python run_release_tests.py --olex2-dir $rd --olex2-exe D:\git\olex2\build\msvc-2026\olex2\x64\Debug\exe\olex2.exe --pythonhome C:\Users\florian\AppData\Local\Python\pythoncore-3.12-64 --salted-model E:\Model_V6 --keep-scratch --no-fetch --timeout 36000 --data-dir D:\devel\olex2-samples --cases \$cases --out-dir '$out' 2>&1 | Select-String -Pattern '^(PASS|FAIL|SKIP)' | ForEach-Object { \$_.Line.Substring(0, [Math]::Min(230, \$_.Line.Length)) }" > "$S/$name.w$k.harness.txt" &
  k=$((k+1))
done
wait
cat "$S"/$name.w*.harness.txt > "$S/$name.harness.txt"
: > "$S/$name.match.txt"
for c in "$S"/.chunk_${name}_*; do
  k=$(( 10#$(echo "$c" | sed "s/.*_//") + 1 ))
  bash "$S/matchall.sh" /d/devel/olex2-release-out/${name}_w$k $(cat "$c") >> "$S/$name.match.txt" 2>&1
done
awk '/^## /{if(s)print s, m, w, n; s=$2; w=0; n=0; m="-"} /matched/{m=$NF} /<-- NOISE/{n++} /<-- [A-Z][a-z]?$/{w++} END{print s, m, w, n}' "$S/$name.match.txt" > "$S/$name.summary.txt"
awk '{split($2,a,"/"); if(a[1]>=0){M+=a[1];N+=a[2]}; W+=$3; Z+=$4} END{print "matched",M"/"N,"mistyped",W,"noise",Z}' "$S/$name.summary.txt"
