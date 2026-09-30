"""The release gate: every NoSpherA2 pathway a major Olex2 update must keep.

One case per pathway on Florian's list (16 Sep 2026), each starting from a
fresh copy of a sample and each returning one line of facts - R1 to four
decimals, the table name, cycle, atom and electron counts, flags - and never
a timing, a path or a version label. run_pipeline.py writes those lines,
sorted by case name, to OLEX2_TEST_RELEASE_LOG; the driver in the
olex2-release-tests repository compares that file with the golden of the
same tier (compare_golden.py: exact line, else integers exact and floats
within 0.0005).

  quick tier (minutes): Thakkar IAM, OCC r2SCAN/def2-SVP, pTB with ECPs on a
      Pd complex, SALTED Model_V6 (model not shipped: SKIP without it),
      discambMATTS and the "Get discambMATTS" entry when the exe is missing,
      four partitionings, XCW with two lambda steps and its GUI block, the
      grown asymmetric unit (water, Mn(II)), properties as cube files, a
      two-domain HKLF 5 twin and the same as HKLF 4 + TWIN card (epoxide,
      synthetic, BASF and -twin), the
      Hybrid mode on ZP2 (molecule 1 as PART 1, molecule 2 as PART 2) with
      discambMATTS and pTB.
  full tier (OLEX2_TEST_FULL=1, needs ORCA): the quick tier plus ORCA on an
      organic molecule, ORCA with ECPs on malbac, embedded (MOL-CRYSTAL-QMMM)
      ORCA on epoxide, Hybrid pTB + ORCA and discambMATTS + ORCA on ZP2.

Nothing here is Tonto. A case whose program or model is not installed raises
SkipTest with the reason; the driver prints those as a banner and exits 2 so
a skip is never mistaken for a pass.

  OLEX2_TEST_FULL          1 registers the full-tier cases as well
  OLEX2_TEST_CASES         comma-separated case names, default all
  OLEX2_TEST_FLINT_SAMPLES  extra flint_<sample> cases beyond the four shipped
                           and the eight carbon-free cod_<id> ones (group_flint.INORGANIC)
  OLEX2_TEST_SALTED_MODEL  directory holding the .salted model (read by
                           group_nsa2_matrix.t_enable)
  OLEX2_TEST_NCPUS, OLEX2_TEST_MEM, OLEX2_TEST_SAMPLE_DIR as for the matrix
"""
from __future__ import absolute_import, division, print_function

import glob
import os
import shutil
import struct
import time

import olx
from olexFunctions import OV

from variableFunctions import nsa2_get_param, nsa2_set_param

from pipeline_tests import (macro, SkipTest, atom_count, clear_r1,
                            r1_of_last_refinement)
import group_nsa2_matrix as matrix
from group_nsa2_matrix import (SAMPLES, _prepare, _set_aspherical,
                               _defaults_for, _aspherical_run, _check_moved)
import group_nosphera2 as base
import group_flint
from group_nosphera2 import (_source_name, _refine_capturing, _reason)

GROUP = "release"

# OCC (method, basis) pairs; one case each. Extend here, the golden grows.
OCC_METHOD_BASIS = [("r2scan", "def2-SVP")]

# partitioning schemes of the NoSpherA2 grid, one case each, all with pTB
PARTITIONS = ["Hirshfeld", "MBIS", "EMBIS", "TFVC"]

# Hybrid mode on ZP2: (part 1 generator, part 2 generator, tier). The label
# is what _source_name resolves; molecule 1 is PART 1, molecule 2 PART 2.
# discambMATTS first and pTB second is the order that lost pTB's ECP core
# electrons before the -mtc_ECP fix (utilities.cuqct_tsc read the part
# settings by position, not by PART number, once a discamb part was taken
# out of the list).
HYBRID_COMBOS = [
  ("DiscaMB", "pTB",  "quick"),
  ("pTB",     "ORCA", "full"),
  ("DiscaMB", "ORCA", "full"),
]

XCW_TIMEOUT = 600.0     # seconds for two lambda steps on epoxide
CUBE_TIMEOUT = 180.0    # seconds for three cubes at 0.5 A


def _env_list(name):
  raw = os.environ.get(name, "")
  return [s.strip() for s in raw.split(",") if s.strip()]


def register(suite):
  full = os.environ.get("OLEX2_TEST_FULL", "") == "1"
  wanted = _env_list("OLEX2_TEST_CASES")
  if [w.lower() for w in wanted] == ["all"]:
    wanted = []
  cases = list(CASES)
  known = [c[0] for c in cases]
  for w in wanted:
    if w not in known:
      suite.run(GROUP, w, _unknown_case, w, known)
  suite.run(GROUP, "setup", t_setup)
  for name, tier, fn, kwds in cases:
    if wanted and name not in wanted:
      continue
    if tier == "full" and not full:
      # not registered at all: the quick golden has no line for it
      continue
    suite.run(GROUP, name, fn, suite, **kwds)


def _unknown_case(name, known):
  raise AssertionError("no such case (have %s)" % ", ".join(known))


def t_setup():
  """What the matrix switches on, plus the XCW block; nothing machine-specific
  in the returned line."""
  base.t_present()
  matrix.t_enable()
  OV.SetParam('user.NoSpherA2.show_XCW', True)
  return ""


# ---------------------------------------------------------------------------
# facts and logs
# ---------------------------------------------------------------------------

def _facts(*pairs):
  return " ".join("%s=%s" % (k, v) for k, v in pairs)


def _r1_facts(tsc, cycles, r_sph, r_asp):
  return [("tsc", os.path.basename(tsc)), ("cycles", cycles),
          ("R1_sph", "%.4f" % r_sph), ("R1_asp", "%.4f" % r_asp)]


def _partition_log(folder):
  """NoSpherA2.log of the last table run, renamed <name>.partitionlog."""
  logs = glob.glob(os.path.join(folder, "*.partitionlog"))
  if not logs:
    raise AssertionError("no .partitionlog beside the model in %s"
                         % os.path.basename(folder))
  logs.sort(key=os.path.getmtime)
  return open(logs[-1], "r", errors="ignore").read()


def _electron_counts(text):
  """protons, electrons, ECP electrons as NoSpherA2 printed them."""
  out = {"protons": None, "electrons": None, "ecp": 0}
  for line in text.splitlines():
    if line.startswith("Number of protons:"):
      out["protons"] = int(float(line.split(":")[1]))
    elif line.startswith("Number of electrons:"):
      out["electrons"] = float(line.split(":")[1])
    elif line.startswith("Number of ECP electrons:"):
      out["ecp"] = int(float(line.split(":")[1]))
  if out["protons"] is None or out["electrons"] is None:
    raise AssertionError("NoSpherA2 did not report its electron counts")
  return out


def _partition_names(text):
  """The columns of 'Table of Charges in electrons' - the scheme(s) used."""
  lines = text.splitlines()
  for i, line in enumerate(lines):
    if line.startswith("Table of Charges in electrons") and i + 1 < len(lines):
      cols = lines[i + 1].split()
      if cols and cols[0] == "Atom":
        return cols[1:]
  return []


def _ecp_facts(counts, charge=0):
  expected = counts["protons"] - counts["ecp"] - charge
  if abs(counts["electrons"] - expected) > 0.01:
    raise AssertionError("protons %d - ECP %d - charge %d = %d, but the wavefunction "
                         "holds %.2f electrons" % (counts["protons"], counts["ecp"],
                                                    charge, expected, counts["electrons"]))
  return [("protons", counts["protons"]), ("electrons", "%.1f" % counts["electrons"]),
          ("ecp", counts["ecp"])]


def _har1(suite, sample, backend, source, setup=None, check=_check_moved):
  """Fresh copy, spherical, one aspherical cycle; the folder and the facts."""
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  if setup:
    setup()
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  check(backend, r_sph, r_asp)
  return folder, _r1_facts(tsc, cycles, r_sph, r_asp)


def _check_not_worse(backend, r_sph, r_asp):
  """For a table that is meant to reproduce the IAM: it may not move R1 much."""
  if r_asp > r_sph + 0.01:
    raise AssertionError("%s: R1 %.4f is worse than spherical %.4f"
                         % (backend, r_asp, r_sph))


class _params(object):
  """Set NoSpherA2 parameters for one case and put them back afterwards."""

  def __init__(self, **kwds):
    self.kwds = kwds
    self.saved = {}

  def __enter__(self):
    for k, v in self.kwds.items():
      self.saved[k] = nsa2_get_param(k)
      nsa2_set_param(k, v)
    return self

  def __exit__(self, *exc):
    for k, v in self.saved.items():
      try:
        nsa2_set_param(k, v)
      except Exception:
        pass
    return False


# ---------------------------------------------------------------------------
# the cases
# ---------------------------------------------------------------------------

def c_iam_thakkar(suite):
  source = _source_name("Thakkar IAM")
  folder, facts = _har1(suite, "epoxide", "Thakkar IAM", source,
                        check=_check_not_worse)
  return _facts(*facts)


def c_occ(suite, method, basis):
  source = matrix._available("OCC")
  methods = os.path.join(OV.BaseDir(), "occ", "share", "methods", "dft_methods.json")
  if os.path.isfile(methods) and method.lower() not in open(methods).read().lower():
    raise SkipTest("occ does not know the method %s (%s)" % (method, methods))

  def setup():
    nsa2_set_param('method', method)
    nsa2_set_param('basis_name', basis)
  folder, facts = _har1(suite, "epoxide", "OCC", source, setup)
  toml = os.path.join(folder, "olex2", "Wfn_job", "epoxide.toml")
  if not os.path.isfile(toml):
    raise AssertionError("no occ input written (%s)" % toml)
  text = open(toml).read().lower()
  for key, value in (("method", method), ("basis", basis)):
    if ('%s = "%s"' % (key, value.lower())) not in text:
      raise AssertionError("occ input does not carry %s = %s" % (key, value))
  return _facts(*(facts + [("method", method), ("basis", basis)]))


def c_ptb_ecp(suite):
  source = _source_name("pTB")
  folder, facts = _har1(suite, "malbac", "pTB", source)
  counts = _electron_counts(_partition_log(folder))
  if counts["ecp"] < 28:
    raise AssertionError("pTB on a Pd complex reports %d ECP electrons, Pd alone has 28"
                         % counts["ecp"])
  return _facts(*(facts + _ecp_facts(counts)))


def c_salted(suite, sample):
  # only the model the caller names counts: a directory remembered in the
  # user settings from an earlier session would test somebody else's model
  model_dir = os.environ.get("OLEX2_TEST_SALTED_MODEL", "").strip()
  if not model_dir:
    raise SkipTest("no SALTED model given - pass --salted-model <dir holding Model_V6>")
  models = sorted(glob.glob(os.path.join(model_dir, "*.salted")))
  if not models:
    raise SkipTest("no .salted file in %s" % model_dir)
  source = matrix._available("SALTED")
  folder, facts = _har1(suite, sample, "SALTED", source)
  return _facts(*(facts + [("model", os.path.basename(models[0]))]))


def c_discamb(suite):
  source = _source_name("DiscaMB")
  folder, facts = _har1(suite, "sucrose", "DiscaMB", source)
  return _facts(*facts)


def c_discamb_gui_get_entry(suite):
  """The source list offers 'Get discambMATTS' exactly once when the
  executable cannot be found, and choosing it opens the DiscaMB page.

  olx.file.Which searches the base directory too, so on a run directory that
  ships discambMATTS2tsc the entry never appears by itself: Which and Shell
  are replaced for the duration of the check, everything is put back after.
  """
  from NoSpherA2.NoSpherA2 import NoSpherA2_instance as nsp2, change_tsc_generator
  opened = []
  real_which, real_shell = olx.file.Which, olx.Shell
  saved_softwares, saved_exe = nsp2.softwares, nsp2.discamb_exe

  def no_which(*args, **kwds):
    return None

  def shell(url, *args, **kwds):
    opened.append(str(url))
    return ""
  olx.file.Which, olx.Shell = no_which, shell
  try:
    nsp2.softwares = ""
    nsp2.setup_discamb()
    entries = [s.strip() for s in nsp2.softwares.split(";") if s.strip()]
    n_get = entries.count("Get discambMATTS")
    if n_get == 0:
      raise AssertionError("no 'Get discambMATTS' entry although the exe is missing "
                           "(list: %s)" % "; ".join(entries) or "empty")
    if any(e.startswith("discamb") for e in entries):
      raise AssertionError("discamb is still offered without an executable: %s"
                           % "; ".join(entries))
    change_tsc_generator("  Get discambMATTS")
    if not opened or "discamb" not in opened[-1].lower():
      raise AssertionError("choosing the entry did not open the DiscaMB page (%s)"
                           % (opened or "nothing"))
  finally:
    olx.file.Which, olx.Shell = real_which, real_shell
    nsp2.softwares, nsp2.discamb_exe = saved_softwares, saved_exe
  return _facts(("get_entries", n_get), ("opened", "discamb"))


def c_orca(suite):
  source = _source_name("ORCA")
  folder, facts = _har1(suite, "epoxide", "ORCA", source)
  return _facts(*facts)


def c_orca_ecp(suite):
  source = _source_name("ORCA")

  def setup():
    nsa2_set_param('basis_name', 'ECP-def2-SVP')
    nsa2_set_param('NoSpherA2_ECP', 1)
  folder, facts = _har1(suite, "malbac", "ORCA", source, setup)
  counts = _electron_counts(_partition_log(folder))
  if counts["ecp"] != 28:
    raise AssertionError("def2 ECP on Pd is 28 core electrons, NoSpherA2 read %d"
                         % counts["ecp"])
  return _facts(*(facts + _ecp_facts(counts) + [("basis", "ECP-def2-SVP")]))


def c_partition(suite, scheme):
  source = _source_name("pTB")
  with _params(NoSpherA2_Partition=scheme):
    folder, facts = _har1(suite, "epoxide", "pTB", source)
    text = _partition_log(folder)
  names = _partition_names(text)
  if scheme not in names:
    raise AssertionError("the charge table is headed %s, not %s - the -%s flag did "
                         "not reach NoSpherA2" % (names or "nothing", scheme, scheme))
  counts = _electron_counts(text)
  if counts["ecp"] <= 0:
    raise AssertionError("pTB ran without -ECP 3: no ECP electrons reported")
  return _facts(*(facts + [("partition", scheme), ("ecp", counts["ecp"])]))


def _gui_html(make):
  """The html a GUI builder returns; olex2c has no xgrid, which the button
  state check in htmlTools asks for, so the console build gets a stand-in
  that answers 'not visible' for the duration of the call."""
  if hasattr(olx, "xgrid"):
    return make()

  class _no_grid(object):
    def Visible(self):
      return "false"
  olx.xgrid = _no_grid()
  try:
    return make()
  finally:
    del olx.xgrid


def c_xcw(suite):
  """Two lambda steps of the X-ray constrained wavefunction on epoxide, the
  GUI block that drives it, and a refinement against the chosen table."""
  import xcw
  from hybrid_GUI import make_XCW_GUI
  folder, r_sph = _prepare(suite, "epoxide")
  charge, mult, doc = SAMPLES["epoxide"]
  nsa2_set_param('charge', str(charge))
  nsa2_set_param('multiplicity', str(mult))
  job = os.path.join(folder, "olex2", "XCW")
  with _params(**{"XCW.start": "0.0", "XCW.step_size": "0.01", "XCW.end": "0.01",
                  "XCW.basis_name": "def2-SVP", "XCW.params": 0}):
    xcw.xcw_run()
    if xcw._state["thread"] is None:
      raise AssertionError("xcw_run did not start (%s)" % (xcw._state["message"] or "no message"))
    t0 = time.time()
    while xcw._state["thread"] is not None and xcw._state["thread"].is_alive():
      if time.time() - t0 > XCW_TIMEOUT:
        xcw.xcw_stop()
        raise AssertionError("XCW did not finish within %.0f s" % XCW_TIMEOUT)
      time.sleep(2.0)
    # the Schedule()d callback never fires while python holds the interpreter
    xcw.xcw_finished()
    if xcw._state["message"]:
      raise AssertionError("XCW: %s" % xcw._state["message"])
    tscbs = xcw.xcw_tscb_list()
    if tscbs != "NA2_0000000.tscb;NA2_0010000.tscb":
      raise AssertionError("expected the tables of lambda 0 and 0.01, got %r" % tscbs)
    status = xcw.xcw_status()
    html = _gui_html(make_XCW_GUI)
    for needle in ("NoSpherA2_XCW_run", "NoSpherA2_XCW_stop", "NoSpherA2_XCW_use",
                   "NA2_0010000.tscb"):
      if needle not in html:
        raise AssertionError("the XCW GUI block lacks %s" % needle)
    xcw.xcw_use_tscb()
    source = str(nsa2_get_param('source'))
    if source != "epoxide_xcw_0010000.tscb":
      raise AssertionError("Use did not make the table the source (source=%r)" % source)
    if not os.path.isfile(os.path.join(folder, source)):
      raise AssertionError("%s was not copied beside the structure" % source)
    _set_aspherical("epoxide", source)
    clear_r1()
    said = _refine_capturing()
    r_xcw = r1_of_last_refinement()
    if r_xcw is None or abs(r_xcw - r_sph) < 1e-5:
      raise AssertionError("refinement against %s left R1 at %s%s"
                           % (source, r_xcw, _reason(said)))
  return _facts(("tscbs", tscbs), ("status", '"%s"' % status), ("gui", "run,stop,use"),
                ("source", source), ("R1_sph", "%.4f" % r_sph), ("R1_xcw", "%.4f" % r_xcw))


def _lattice_atoms():
  """Atoms in the lattice, which is what grow adds to; the asymmetric unit
  count does not move."""
  n = 0
  for i in range(int(olx.xf.latt.GetFragmentCount())):
    names = str(olx.xf.latt.GetFragmentAtoms(i))
    n += len([s for s in names.split(",") if s.strip()])
  return n


def c_grown(suite):
  """The asymmetric unit grown to the whole molecule before the table."""
  source = _source_name("pTB")
  folder, r_sph = _prepare(suite, "water")
  # the refinement left Q peaks; NoSpherA2 kills them before it writes the
  # xyz, and olex2c's kill rebuilds the lattice from the asymmetric unit,
  # which would silently undo the grow. Killing them first keeps the two
  # builds on the same path (the GUI's kill leaves a grown lattice alone).
  macro("kill $Q")
  n_asu = _lattice_atoms()
  macro("grow")
  n_grown = _lattice_atoms()
  if n_grown <= n_asu:
    raise AssertionError("grow left the atom count at %d" % n_grown)
  if str(olx.xf.latt.IsGrown()).lower() != "true":
    raise AssertionError("the lattice does not report itself as grown")
  _set_aspherical("water", source)
  cycles, tsc, said = _aspherical_run("water", "pTB", source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved("pTB", r_sph, r_asp)
  # whether the lattice is still grown after the refinement differs between
  # the GUI and olex2c, so it is not part of the golden line
  if str(olx.xf.latt.IsGrown()).lower() == "true":
    macro("fuse")
  macro("kill $Q")
  n_fused = _lattice_atoms()
  if n_fused != n_asu:
    raise AssertionError("fuse gives %d atoms, the asymmetric unit had %d" % (n_fused, n_asu))
  return _facts(*(_r1_facts(tsc, cycles, r_sph, r_asp)
                  + [("atoms_asu", n_asu), ("atoms_grown", n_grown),
                     ("multiplicity", SAMPLES["water"][1])]))


HKLF5_LAW = ((1, 0, 0), (0, -1, 0), (0, 0, -1))   # a twofold about a: pseudo-merohedral in P2(1)/c


def _write_hklf5(hkl, law, fraction=0.2, batches=True):
  """Rewrite an HKLF 4 file (monoclinic, merged) as two-component HKLF 5: every reflection whose
  image under the law was measured becomes an overlapped pair (its image on a
  batch -2 line, itself on the batch 1 line that closes the group) with the
  image's intensity mixed in. Synthetic, so it tests the plumbing, not the
  physics. Returns how many reflections were paired. batches=False keeps the
  file HKLF 4 with the same contamination, the case for a TWIN card."""
  rows = []
  for line in open(hkl):
    if len(line) < 28:
      break
    h, k, l = int(line[0:4]), int(line[4:8]), int(line[8:12])
    if (h, k, l) == (0, 0, 0):
      break
    rows.append(((h, k, l), float(line[12:20]), float(line[20:28])))
  # merged data: the image is stored under one of its 2/m equivalents
  data = {}
  for (h, k, l), i, s in rows:
    for eq in ((h, k, l), (-h, k, -l), (-h, -k, -l), (h, -k, l)):
      data.setdefault(eq, (i, s))
  out, paired = [], 0
  for hkl_, i, s in rows:
    image = tuple(sum(law[r][c] * hkl_[c] for c in range(3)) for r in range(3))
    if image in data and image != hkl_:
      # the group's one observation is read from the closing line
      i2, s2 = data[image]
      if batches:
        out.append("%4d%4d%4d%8.2f%8.2f%4d" % (image + (0.0, 0.0, -2)))
      i, s = (1 - fraction) * i + fraction * i2, (s * s + (fraction * s2) ** 2) ** 0.5
      paired += 1
    out.append("%4d%4d%4d%8.2f%8.2f" % (hkl_ + (i, s)) + ("%4d" % 1 if batches else ""))
  out.append("%4d%4d%4d%8.2f%8.2f" % (0, 0, 0, 0, 0) + ("%4d" % 0 if batches else ""))
  open(hkl, "w").write("\n".join(out) + "\n")
  return paired


def c_hklf5(suite):
  """Non-merohedral twin data (HKLF 5 with batch numbers) through the whole
  route: BASF refines, the law read off the file reaches NoSpherA2 as -twin,
  the table covers the second domain, the aspherical cycle runs."""
  source = _source_name("pTB")
  folder = matrix.sample_copy(suite, "epoxide")
  hkl = [os.path.join(folder, f) for f in os.listdir(folder) if f.lower().endswith(".hkl")][0]
  paired = _write_hklf5(hkl, HKLF5_LAW)
  if paired < 100:
    raise AssertionError("only %d reflections have a measured twin image" % paired)
  matrix._load_model(folder, matrix._model_file(folder))
  macro("AddIns HKLF 5")
  macro("AddIns BASF 0.2")
  macro("file")
  OV.SetParam('snum.refinement.update_weight', False)
  clear_r1()
  base._refine_spherically()
  r_sph = r1_of_last_refinement()
  model = OV.GetRefinementModel(False)
  if model['hklf']['value'] != 5:
    raise AssertionError("the model reports HKLF %s after the spherical refinement" % model['hklf']['value'])
  basf = [float(b) for b in model['hklf'].get('basf', [])]
  if not basf or not 0.0 < basf[0] < 1.0:
    raise AssertionError("BASF did not refine to a fraction: %r" % basf)
  _set_aspherical("epoxide", source)
  cycles, tsc, said = _aspherical_run("epoxide", "pTB", source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved("pTB", r_sph, r_asp)
  from cctbx_olex_adapter import get_table_fallback_atoms
  if get_table_fallback_atoms():
    raise AssertionError("atoms refined spherically beside the table: %s" % get_table_fallback_atoms())
  return _facts(*(_r1_facts(tsc, cycles, r_sph, r_asp)
                  + [("paired", paired), ("basf", "%.3f" % basf[0])]))


def c_twin_hklf4(suite):
  """Merohedral twin as SHELX writes it: HKLF 4 data with a TWIN card and BASF.
  The law (and its powers) reach NoSpherA2 as -twin, so the table covers R h
  and the aspherical cycle runs."""
  source = _source_name("pTB")
  folder = matrix.sample_copy(suite, "epoxide")
  hkl = [os.path.join(folder, f) for f in os.listdir(folder) if f.lower().endswith(".hkl")][0]
  paired = _write_hklf5(hkl, HKLF5_LAW, batches=False)
  if paired < 100:
    raise AssertionError("only %d reflections have a measured twin image" % paired)
  matrix._load_model(folder, matrix._model_file(folder))
  macro("AddIns TWIN %s 2" % " ".join(str(x) for row in HKLF5_LAW for x in row))
  macro("AddIns BASF 0.2")
  macro("file")
  OV.SetParam('snum.refinement.update_weight', False)
  clear_r1()
  base._refine_spherically()
  r_sph = r1_of_last_refinement()
  model = OV.GetRefinementModel(False)
  if model['hklf']['value'] != 4 or 'twin' not in model:
    raise AssertionError("expected HKLF 4 with a TWIN card, got HKLF %s, twin %r"
                         % (model['hklf']['value'], model.get('twin')))
  basf = [float(b) for b in model['twin'].get('basf', [])]
  if not basf or not 0.0 < basf[0] < 1.0:
    raise AssertionError("BASF did not refine to a fraction: %r" % basf)
  _set_aspherical("epoxide", source)
  cycles, tsc, said = _aspherical_run("epoxide", "pTB", source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved("pTB", r_sph, r_asp)
  from cctbx_olex_adapter import get_table_fallback_atoms
  if get_table_fallback_atoms():
    raise AssertionError("atoms refined spherically beside the table: %s" % get_table_fallback_atoms())
  return _facts(*(_r1_facts(tsc, cycles, r_sph, r_asp)
                  + [("paired", paired), ("basf", "%.3f" % basf[0])]))


def c_orca_qmmm(suite):
  """ORCA inside the crystal field (MOL-CRYSTAL-QMMM), the embedded route."""
  source = _source_name("ORCA")
  if source not in ("ORCA 5.0", "ORCA 6.0", "ORCA 6.1"):
    raise SkipTest("embedded ORCA needs the ORCA 5.0/6.0/6.1 label, this one is %r" % source)
  with _params(ORCA_USE_CRYSTAL_QMMM=True, ORCA_CRYSTAL_QMMM_TYPE="Mol",
               ORCA_CRYSTAL_QMMM_RADIUS="6.0"):
    folder, facts = _har1(suite, "epoxide", "ORCA", source)
  inp = os.path.join(folder, "olex2", "Wfn_job", "epoxide.inp")
  if not os.path.isfile(inp):
    raise AssertionError("no ORCA input written (%s)" % inp)
  text = open(inp, errors="ignore").read()
  for needle in ("MOL-CRYSTAL-QMMM", "%qmmm"):
    if needle not in text:
      raise AssertionError("the ORCA input has no %s block - not the embedded route" % needle)
  return _facts(*(facts + [("embedding", "MOL-CRYSTAL-QMMM"), ("radius", "6.0")]))


def _fragment_labels():
  """The atom labels of every fragment of the lattice, one list each."""
  out = []
  for i in range(int(olx.xf.latt.GetFragmentCount())):
    names = str(olx.xf.latt.GetFragmentAtoms(i))
    out.append([s.strip() for s in names.split(",") if s.strip()])
  return out


def _assign_parts(anchor):
  """The fragment holding the anchor atom becomes PART 1, everything else
  PART 2, the way a user would do it with two selections and the part
  button. Returns the atom count of each part."""
  macro("kill $Q")
  frags = _fragment_labels()
  first = [f for f in frags if anchor in f]
  if len(first) != 1:
    raise AssertionError("%s is in %d fragments, expected one (fragments: %s)"
                         % (anchor, len(first), [len(f) for f in frags]))
  rest = [a for f in frags if anchor not in f for a in f]
  if not rest:
    raise AssertionError("only one fragment, nothing left for PART 2")
  macro("part 1 " + " ".join(first[0]))
  macro("part 2 " + " ".join(rest))
  parts = sorted(OV.ListParts() or [])
  if parts != [1, 2]:
    raise AssertionError("the model carries parts %s after the assignment, not [1, 2]" % parts)
  counts = {1: 0, 2: 0}
  for i in range(int(olx.xf.au.GetAtomCount())):
    p = int(olx.xf.au.GetAtomPart(i))
    if p in counts:
      counts[p] += 1
  return counts


def _ecp_total(text):
  """Sum of every 'Number of ECP electrons' line - one per wavefunction in
  a -mtc run."""
  n = 0
  for line in text.splitlines():
    if line.startswith("Number of ECP electrons:"):
      n += int(float(line.split(":")[1]))
  return n


def c_hybrid(suite, part1, part2):
  """Hybrid mode: one generator per PART, the tables partitioned or merged.

  The part settings are written the way the Hybrid GUI writes them - the
  software combos are built from the source list, so their values carry
  its two leading spaces, and discambMATTS is the exe name. That spelling
  is the point: the run must take the GUI's values, not a tidied copy.
  """
  labels = {"DiscaMB": _source_name("DiscaMB"), part1: _source_name(part1),
            part2: _source_name(part2)}
  discamb = labels["DiscaMB"]
  folder, r_sph = _prepare(suite, "ZP2")
  counts = _assign_parts("N11")
  settings = {}
  for n, gen in ((1, part1), (2, part2)):
    settings["Hybrid.software_Part%d" % n] = "  " + labels[gen]
    settings["Hybrid.basis_name_Part%d" % n] = "def2-SVP"
    settings["Hybrid.method_Part%d" % n] = "B3LYP"
    settings["Hybrid.charge_Part%d" % n] = "0"
    settings["Hybrid.multiplicity_Part%d" % n] = "1"
    settings["Hybrid.Relativistic_Part%d" % n] = False
    settings["Hybrid.ORCA_SCF_Conv_Part%d" % n] = "NoSpherA2SCF"
    settings["Hybrid.ORCA_SCF_Strategy_Part%d" % n] = "NormalConv"
  _set_aspherical("ZP2", "Hybrid")
  with _params(**settings):
    cycles, tsc, said = _aspherical_run("ZP2", "Hybrid", "Hybrid", folder)
  r_asp = r1_of_last_refinement()
  _check_moved("Hybrid", r_sph, r_asp)
  # the table the refinement used, not the newest file: a merge leaves the
  # per-part tables beside the total
  used = str(nsa2_get_param('file') or "")
  if not used:
    raise AssertionError("no tsc file recorded after the hybrid run")
  gens = [g for g in (part1, part2) if labels[g] != discamb]
  facts = _r1_facts(used, cycles, r_sph, r_asp)
  facts += [("parts", "1,2"), ("atoms_part1", counts[1]), ("atoms_part2", counts[2]),
            ("part1", part1), ("part2", part2)]
  # every wavefunction part is a Part_<n> job with its output copied beside
  # the model as <name>_part<n>.<ext>; every discamb part left <name>_part_<n>.tsc
  wfn_copies = sorted(os.path.basename(p) for p in glob.glob(os.path.join(folder, "ZP2_part[12].*")))
  part_tables = sorted(os.path.basename(p) for p in glob.glob(os.path.join(folder, "ZP2_part_[12].tsc*")))
  n_discamb = len([g for g in (part1, part2) if labels[g] == discamb])
  if len(wfn_copies) != len(gens):
    raise AssertionError("%d wavefunction parts expected, found %s" % (len(gens), wfn_copies))
  if len(part_tables) != n_discamb:
    raise AssertionError("%d discamb part tables expected, found %s" % (n_discamb, part_tables))
  if n_discamb and "_total.tsc" not in os.path.basename(used):
    raise AssertionError("a discamb part must end in a merged _total.tsc(b), the refinement used %s" % used)
  facts += [("wfn", ";".join(wfn_copies) or "-"), ("tables", ";".join(part_tables) or "-")]
  if gens:
    ecp = _ecp_total(_partition_log(folder))
    if "pTB" in gens and ecp <= 0:
      raise AssertionError("the pTB part ran without its ECP core electrons (-mtc_ECP)")
    facts.append(("ecp", ecp))
  return _facts(*facts)


def _cube_dims(path):
  # NoSpherA2 writes binary .cubeb since 21 Sep 2026: magic, version, offset, na, size[3]
  with open(path, "rb") as f:
    if f.read(8) == b"NSA2CUBE":
      f.read(8)
      na, nx, ny, nz = struct.unpack("<4i", f.read(16))
      return abs(na), [nx, ny, nz]
  lines = open(path, errors="ignore").read().splitlines()
  if len(lines) < 7:
    raise AssertionError("%s is not a cube file" % os.path.basename(path))
  natoms = abs(int(lines[2].split()[0]))
  dims = [int(lines[i].split()[0]) for i in (3, 4, 5)]
  return natoms, dims


def _settle(cubes, timeout=CUBE_TIMEOUT):
  """Wait for a properties run to finish: calculate_cubes does not wait, so a
  cube is done only once it is on disk and has stopped growing. Hands back the
  paths with the .cubeb NoSpherA2 writes resolved."""
  log = os.path.join(os.path.dirname(cubes[0]), "NoSpherA2_cube.log")
  t0 = time.time()
  sizes = {}
  while True:
    done = True
    cubes = [c + "b" if os.path.isfile(c + "b") else c for c in cubes]
    for c in cubes:
      size = os.path.getsize(c) if os.path.isfile(c) else -1
      if size <= 0 or sizes.get(c) != size:
        done = False
      sizes[c] = size
    if done:
      return cubes
    if time.time() - t0 > timeout:
      tail = ""
      if os.path.isfile(log):
        tail = " - log: " + " | ".join(open(log, errors="ignore").read().splitlines()[-3:])
      missing = [os.path.basename(c) for c in cubes if not os.path.isfile(c)]
      raise AssertionError("cubes not written within %.0f s (missing %s)%s"
                           % (timeout, ", ".join(missing) or "none, still growing", tail))
    time.sleep(2.0)


def c_cubes(suite):
  """Properties from the pTB wavefunction: Laplacian, ELI-D and the static
  deformation density as cube files, offered afterwards as map types."""
  import cubes_maps
  source = _source_name("pTB")
  folder, facts = _har1(suite, "epoxide", "pTB", source)
  name = OV.ModelSrc()
  wfn = os.path.join(folder, name + ".xtb")
  if not os.path.isfile(wfn):
    inner = os.path.join(folder, "olex2", "Wfn_job", name + ".xtb")
    if not os.path.isfile(inner):
      raise AssertionError("pTB left no %s.xtb beside the model or in the job folder" % name)
    shutil.copy(inner, wfn)
  wanted = ["lap", "eli", "def"]
  # products live where cubes_maps puts them (olex2/NoSpherA2/<key>/<wfn key>/ since 21 Sep 2026), as .cubeb
  cubes = [cubes_maps._p("%s_%s.cube" % (name, w)) for w in wanted]
  for c in cubes:
    for f in (c, c + "b"):
      if os.path.isfile(f):
        os.remove(f)
  os.chdir(folder)
  with _params(**{"Property_Lap": True, "Property_Eli": True, "Property_DEF": True,
                  "Property_Elf": False, "Property_RDG": False, "Property_ESP": False,
                  "Property_MO": False, "Property_ATOM": False, "Property_all_MOs": False,
                  "map.radius": "1.0", "map.resolution": "0.5"}):
    cubes_maps.calculate_cubes()
  cubes = _settle(cubes)
  natoms, dims = _cube_dims(cubes[0])
  for c in cubes[1:]:
    n2, d2 = _cube_dims(c)
    if (n2, d2) != (natoms, dims):
      raise AssertionError("%s has a different grid than %s"
                           % (os.path.basename(c), os.path.basename(cubes[0])))
  types = str(cubes_maps.get_map_types())
  for label in ("Laplacian;", "ELI-D;", "Stat. Def.;"):
    if label not in types:
      raise AssertionError("map type %r not offered although its cube exists (%s)"
                           % (label, types))
  return _facts(*(facts + [("cubes", ",".join(wanted)), ("cube_atoms", natoms),
                           ("grid", "x".join(str(d) for d in dims))]))


# name, tier, function, keyword arguments - the order is the golden's order
CASES = [
  ("iam_thakkar_epoxide",    "quick", c_iam_thakkar, {}),
  ("ptb_ecp_malbac",         "quick", c_ptb_ecp, {}),
  ("salted_v6_sucrose",      "quick", c_salted, {"sample": "sucrose"}),
  ("salted_v6_epoxide",      "quick", c_salted, {"sample": "epoxide"}),
  ("discamb_sucrose",        "quick", c_discamb, {}),
  ("discamb_gui_get_entry",  "quick", c_discamb_gui_get_entry, {}),
  ("xcw_epoxide",            "quick", c_xcw, {}),
  ("grown_water_ptb",        "quick", c_grown, {}),
  ("cubes_epoxide_ptb",      "quick", c_cubes, {}),
  ("hklf5_epoxide_ptb",      "quick", c_hklf5, {}),
  ("twin_hklf4_epoxide_ptb", "quick", c_twin_hklf4, {}),
  ("flint_sucrose",      "quick", group_flint.t_flint, {"sample": "sucrose"}),
  ("flint_epoxide",      "quick", group_flint.t_flint, {"sample": "epoxide"}),
  ("flint_water",        "quick", group_flint.t_flint, {"sample": "water"}),
  ("flint_malbac",       "quick", group_flint.t_flint, {"sample": "malbac"}),
  ("flint_cod_1560875",  "quick", group_flint.t_flint, {"sample": "cod_1560875"}),
  ("flint_cod_2241658",  "quick", group_flint.t_flint, {"sample": "cod_2241658"}),
  ("flint_cod_2108240",  "quick", group_flint.t_flint, {"sample": "cod_2108240"}),
  ("flint_cod_2013004",  "quick", group_flint.t_flint, {"sample": "cod_2013004"}),
  ("flint_cod_2104335",  "quick", group_flint.t_flint, {"sample": "cod_2104335"}),
  ("flint_cod_2208447",  "quick", group_flint.t_flint, {"sample": "cod_2208447"}),
  ("flint_cod_2229150",  "quick", group_flint.t_flint, {"sample": "cod_2229150"}),
  ("flint_cod_2108989",  "quick", group_flint.t_flint, {"sample": "cod_2108989"}),
  ("orca_epoxide",           "full",  c_orca, {}),
  ("orca_ecp_malbac",        "full",  c_orca_ecp, {}),
  ("orca_qmmm_epoxide",      "full",  c_orca_qmmm, {}),
]
for _s in _env_list("OLEX2_TEST_FLINT_SAMPLES"):
  CASES.append(("flint_%s" % _s.lower(), "quick", group_flint.t_flint,
                {"sample": _s}))
for _scheme in PARTITIONS:
  CASES.append(("part_%s_epoxide" % _scheme.lower(), "quick", c_partition,
                {"scheme": _scheme}))
for _p1, _p2, _tier in HYBRID_COMBOS:
  CASES.append(("hybrid_zp2_%s_%s" % (_p1.lower(), _p2.lower()), _tier, c_hybrid,
                {"part1": _p1, "part2": _p2}))
for _method, _basis in OCC_METHOD_BASIS:
  CASES.append(("occ_%s_%s_epoxide" % (_method.lower(), _basis), "quick", c_occ,
                {"method": _method, "basis": _basis}))
# pictures of the 3D view, each from a camera matr fixed
import render_tests
CASES += render_tests.CASES
# the GUI layer: the panels first, on a fresh window, and the log check last
# so it covers everything the refinements printed as well
import gui_tests
CASES = ([(n, "quick", f, {}) for n, f in gui_tests.CASES if n != "gui_log_clean"]
         + CASES
         + [(n, "quick", f, {}) for n, f in gui_tests.CASES if n == "gui_log_clean"])
