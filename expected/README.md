# Golden logs

`release_tests.quick.good` and `release_tests.full.good` are the reference
logs, one line per case (sorted) plus `END <pass> <fail> <skip>`. They were
produced on FLOWOFFICE (Windows 11) on 16 Sep 2026 with the NoSpherA2 bundle of
commit 0442ebf1 (XCW_Test), ORCA 6.1, and SALTED Model_V6.

Re-recorded 30 Sep 2026 (quick 36, full 45): the render cases added, FLINT
after the rename and the whole-COD fixes. Olex2 was a Debug build of 11a9a51b
plus the then-uncommitted `AddObject -t` (surface transparency) - a build
without it fails `render_hirshfeld`/`render_esp_surface`; NoSpherA2 the
rundir-py3 exe of 21 Sep. `flint_cod_2108240` swaps one atom between Al and
Si from run to run (`types=`), so its text may differ while `typed` holds.

Accept a new golden only after the diff is understood:

    python run_release_tests.py ... --keep-scratch      # look at the case
    python run_release_tests.py ... --update-golden     # then accept

A per-platform file `release_tests.<tier>.<win|linux|mac>.good` overrides the
shared one when present; add one only if a platform drifts beyond 0.0005 in R1.

# Reference pictures

`gui/<tier>.<platform>_<W>x<H>/` holds the screenshots of the GUI and render
cases of one tier on one screen size. `compare_gui.py` fails the run when
a picture changed (any pixel off by more than 24 in a channel outside
`mask.txt`), is new or is missing, and paints the differing pixels magenta in
`<out>/gui/diff/`. The render cases' GL pictures (`view_*`, `hirshfeld_*`,
...) are compared and listed as `RENDER` but do not fail: `pict` draws the
whole scene at one of two scales from run to run (ratio about the canvas
aspect, cause not found yet; `<state>.zoom.txt` holds the zoom numbers it
used, unchanged between the two). On a screen size with no folder here the pictures are only
noted; `--update-golden` on the whole tier records them. `mask.txt` holds the
rectangles that print the out dir (title bar, panel header, status bar); the
default was measured on FLOWOFFICE, check it on the first diff of a new screen.
