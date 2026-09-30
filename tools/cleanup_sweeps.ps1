# Cleanup of finished FLINT sweep output. Keeps every log/harness file, removes the bulk
# (per-worker sample staging and Olex2 scratch). Run on the machine named in each block.
# Passes still running are excluded by name; -WhatIf first if unsure.

# ---- SPYRO (idle, all results pulled to FLOWOFFICE) ----
# 1) everything that ran from D:\ before the E:-only rule: whole tree goes
Remove-Item -Recurse -Force D:\devel\olex2-release-out
# 2) E:\ passes: keep the tiny logs, drop samples + scratch
Get-ChildItem E:\olex2-release-out -Directory | ForEach-Object {
  foreach ($d in 'samples', 'scratch') { $p = Join-Path $_.FullName $d; if (Test-Path $p) { Remove-Item -Recurse -Force $p } }
}
# 3) drained job queue
Remove-Item -Force E:\olex2-jobs\done\*.cmd, E:\olex2-jobs\done\*.log -ErrorAction SilentlyContinue

# ---- FLOWOFFICE (fullA/fullB running, extL next, enorm kept: excluded) ----
$keep = 'fullA_w*', 'fullB_w*', 'extL*', 'enorm*'
Get-ChildItem D:\devel\olex2-release-out -Directory |
  Where-Object { $n = $_.Name; -not ($keep | Where-Object { $n -like $_ }) } |
  ForEach-Object {
    foreach ($d in 'samples', 'scratch') { $p = Join-Path $_.FullName $d; if (Test-Path $p) { Remove-Item -Recurse -Force $p } }
  }
