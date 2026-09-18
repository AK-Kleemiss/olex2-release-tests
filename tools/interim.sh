#!/bin/bash
# interim.sh <pass> <nworkers> : match every finished scratch case (all but the newest per worker) while the run is going
S="$(dirname "$0")"; name=$1; nw=$2
: > "$S/$name.interim.match.txt"
for k in $(seq 1 $nw); do
  d=/d/devel/olex2-release-out/${name}_w$k
  done_=$(ls -t "$d/scratch" 2>/dev/null | tail -n +2)
  [ -n "$done_" ] && bash "$S/matchall.sh" "$d" $done_ >> "$S/$name.interim.match.txt" 2>&1
done
awk '/^## /{if(s)print s, m, w, n; s=$2; w=0; n=0; m="-"} /matched/{m=$NF} /<-- NOISE/{n++} /<-- [A-Z][a-z]?$/{w++} END{print s, m, w, n}' "$S/$name.interim.match.txt" > "$S/$name.interim.summary.txt"
awk '{split($2,a,"/"); if(a[1]>=0){M+=a[1];N+=a[2]}; W+=$3; Z+=$4; C++} END{print C, "cases matched",M"/"N,"mistyped",W,"noise",Z}' "$S/$name.interim.summary.txt"
