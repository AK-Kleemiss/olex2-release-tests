# olex2-release-tests

Release gate for the NoSpherA2 paths of Olex2: a fixed list of feature paths is
run inside the *shipped* Olex2 bundle, the result is written as a plain-text log,
and the log is compared with a golden log kept in this repository. Meant to be
run before a major Olex2 release, on Windows, Linux and macOS.

The test code (`olex2_pipeline_tests/`, group `release`) runs inside Olex2
but is not shipped with it: the same files live in the Olex2 SVN under
`util/pyUtil/regression/olex2_pipeline_tests/` without the `olex-install` /
`olex-update` properties, so users never see them. Edit them in the SVN and
copy them here (or the other way round; `diff -r` between the two says which
side is behind). This repository is the whole test: driver, test modules,
comparer, sample structures and goldens.

## What is tested

| case | path exercised | tier |
|---|---|---|
| `iam_thakkar_epoxide` | Thakkar IAM form factors (smallest test) | quick |
| `occ_r2scan_def2-SVP_epoxide` | OCC wavefunction (r2SCAN / def2-SVP) | quick |
| `ptb_ecp_malbac` | pTB from the bundle, ECP electrons counted (Pd) | quick |
| `salted_v6_sucrose`, `salted_v6_epoxide` | SALTED prediction with a Model_V6 (`--salted-model`) | quick |
| `discamb_sucrose` | discambMATTS2tsc from the bundle | quick |
| `discamb_gui_get_entry` | the GUI offers "Get discambMATTS" when the exe is missing | quick |
| `part_hirshfeld_epoxide`, `part_mbis_epoxide`, `part_embis_epoxide`, `part_tfvc_epoxide` | the four partitionings | quick |
| `xcw_epoxide` | XCW with two lambda steps, its GUI, use of the resulting tscb | quick |
| `grown_water_ptb` | refinement on a grown asymmetric unit (Z' = 0.5, multiplicity 6) | quick |
| `cubes_epoxide_ptb` | property cubes (Laplacian, ELI, deformation density) | quick |
| `hybrid_zp2_discamb_ptb` | Hybrid mode, PART 1 = discambMATTS, PART 2 = pTB, tables merged | quick |
| `orca_epoxide` | ORCA B3LYP/def2-SVP | full |
| `orca_ecp_malbac` | ORCA with an ECP basis, ECP electrons counted | full |
| `orca_qmmm_epoxide` | embedded (MOL-CRYSTAL-QMMM) ORCA | full |
| `hybrid_zp2_ptb_orca` | Hybrid mode, PART 1 = pTB, PART 2 = ORCA | full |
| `hybrid_zp2_discamb_orca` | Hybrid mode, PART 1 = discambMATTS, PART 2 = ORCA | full |
| `gui_tabs_and_panels` | every main tab and tool-bar panel renders without an error line | quick |
| `gui_refine_nosphera2` | the refine panel with NoSpherA2 on: its controls and three blocks | quick |
| `gui_sources` | each entry of the source combo rebuilds the NoSpherA2 block | quick |
| `gui_log_clean` | the error lines of the whole Olex2 log, as known kinds (goal: 0) | quick |

Every case runs one HAR cycle and records the spherical and the aspherical R1,
the table used, and the integers that prove the path was really taken (ECP
electrons, atom counts, part numbers, cube grid, ...). Tonto is not part of it.

## Running

```
python run_release_tests.py --olex2-dir /path/to/unpacked/olex2 \
                            --salted-model /path/to/Model_V6 [--full]
```

- Before anything starts the driver runs `git fetch` and fast-forwards this
  clone when GitHub is ahead, so the goldens are the current ones; a clone
  with local changes is not touched but gets a banner. `--no-fetch` skips it.
- `--olex2-dir` is the unpacked bundle (the directory with `olex2.tag`, the
  NoSpherA2 executable and `util/pyUtil`). It is not modified beyond a
  `runonce.release_tests.txm` file that is removed when the run is over.
- Everything else (data directory, refinement scratch, a copy of the samples,
  the log) goes to `--out-dir` (default `release_out/` next to the script).
- The `gui_*` cases (GUI build only, skipped with olex2c) also leave
  `<out-dir>/gui/index.html`: a screenshot of the Olex2 window for every tab,
  panel and source next to the html tree it was rendered from. The golden
  holds what a script can decide; the pictures are for a human look at
  clipped labels, overlapping controls or an empty panel.
- `--full` adds the ORCA cases; ORCA must be on `PATH` or configured in Olex2.
- Linux without a display: `xvfb-run -a python run_release_tests.py ...`.
- Windows developers can use the console build instead of the GUI:
  `--olex2c path/to/olex2c.exe --pythonhome <python the bundle was built for>`.
- `--cases a,b` runs a subset and compares only those lines (reported as a
  subset, not a release result; `--update-golden` refuses a subset).

Exit codes: `0` every case matches the golden; `2` everything that ran matches
but cases were SKIPPED (a banner names the missing program or model); `1` a
case failed, a number drifted, a case is missing or Olex2 did not finish.

## The golden

`expected/release_tests.<tier>.good` is one line per case, sorted by name, and
a final `END <pass> <fail> <skip>`. `compare_golden.py` demands an exact line
or, failing that, the same skeleton with integers equal and floats within
`0.0005` (the R1 tolerance). A drift beyond that is a real change: look at the
case's scratch (`--keep-scratch`) before accepting it with `--update-golden`.
A platform-specific `release_tests.<tier>.<win|linux|mac>.good` is used when
present; none is needed so far.

## Samples

`sample_data/` holds small structures with reflections: `epoxide` (NoSpherA2
test set), `sucrose`, `water`, `ZP2` (Olex2 `sample_data`), `malbac` (Olex2
`sample_data` structure, reflections from the NoSpherA2 test set). See
`NOTICE`.
