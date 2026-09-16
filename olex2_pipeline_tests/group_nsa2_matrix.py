"""NoSpherA2 refinement paths x backends x sample structures.

group_nosphera2 asks one question per backend: does a single aspherical cycle
on sucrose move R1. This group asks the questions a user's session asks after
that one - does the HAR loop iterate and stop, do anisotropic or free hydrogens
survive the loop, can the table of a finished run be refined against again,
does a table-only run leave a table and nothing else, does the fcf-only route
write an fcf, does a second method/basis of the same backend go through, and
does a disordered structure still get one table through `-mtc`.

Every case is `sample x backend x path`, and every case starts from a fresh
copy of the sample, so a path never inherits the previous path's table or
refined model. Paths that need an earlier result (tsc_reuse) produce it
themselves inside the case.

Selection is by environment variable, as for the rest of the pipeline tests:

  OLEX2_TEST_SAMPLES    comma-separated sample names, default "sucrose".
                        Known: epoxide (fast), sucrose, 183, THPP (disorder).
  OLEX2_TEST_BACKENDS   comma-separated backend labels, default all of BACKENDS
  OLEX2_TEST_PATHS      comma-separated path names, default "spherical,har1"
                        "all" runs every path
  OLEX2_TEST_FULL       1 to include the quantum backends and the slow paths
  OLEX2_TEST_NCPUS      cores per job (snum.NoSpherA2.ncpus), default 8 or fewer
  OLEX2_TEST_MEM        GB per job (snum.NoSpherA2.mem), default 4
  OLEX2_TEST_SALTED_MODEL
                        directory holding model.salted; installed into
                        user.NoSpherA2.salted_models_list for the run
  OLEX2_TEST_SAMPLE_DIR extra directory to look for samples in

The samples are looked up in <basedir>/sample_data, <basedir>/samples,
user.sample_dir and OLEX2_TEST_SAMPLE_DIR, in that order: olex2c never runs
the GUI start-up that copies sample_data into the user's samples folder, so a
fresh run directory has only sample_data.

Backends that are behind a switch in the GUI are switched on here for the run:
SALTED needs olex2.debug, OCC needs user.NoSpherA2.show_OCC *at start-up*
(setup_occ_executables is what adds it to the list, so it is called again).
The switches are session-only settings, nothing is written to the user's phil.
"""
from __future__ import absolute_import, division, print_function

import os
import shutil
import sys
import time

import olx
import olex
from olexFunctions import OV

from variableFunctions import nsa2_get_param, nsa2_set_param

from pipeline_tests import (macro, SkipTest, load, model_in, has_hkl,
                            atom_count,
                            clear_r1, r1_of_last_refinement, _r1_raw,
                            Result)
import group_nosphera2 as base
from group_nosphera2 import (_source_name, _offered, _har_log_cycles,
                             _tsc_beside, _refine_capturing, _reason,
                             _refine_spherically, _all_tables,
                             _per_part_tables, _wavefunctions)

# name -> (charge, multiplicity, why it is here)
SAMPLES = {
  "epoxide": (0, 1, "8 non-H atoms, 2082 reflections: the fast case"),
  "sucrose": (0, 1, "the standard HAR case, 45 atoms"),
  "183":     (0, 1, "C/H/F, 36 atoms"),
  "THPP":    (0, 1, "carries PART, the -mtc case"),
  "malbac":  (0, 1, "Pd complex, one PART -1 group: the ECP case"),
  "water":   (0, 6, "Mn(II) on an inversion centre, Z' = 0.5: grown before the table"),
  "ZP2":     (0, 1, "two independent C8H11F2N3O molecules, Z' = 2: the Hybrid/PART case"),
}
DEFAULT_SAMPLES = "sucrose"

# label -> quantum backend (minutes), so behind OLEX2_TEST_FULL
BACKENDS = [
  ("DiscaMB", False),
  ("SALTED",  False),
  ("xTB",     False),
  ("pTB",     False),
  ("ORCA",    True),
  ("OCC",     True),
]
WAVEFUNCTION_BACKENDS = ["xTB", "pTB", "ORCA", "OCC"]

# backend -> a second (method, basis) to run under alt_method. The default is
# B3LYP/def2-SVP where the backend takes a functional at all.
ALT_METHOD = {
  "ORCA":  ("PBE", "def2-SVP"),
  "OCC":   ("PBE", "def2-SVP"),
  "xTB":   ("GFN1", None),
}

# name -> (slow, what it exercises)
PATHS = [
  ("spherical",  False, "olex2.refine baseline without a table"),
  ("har1",       False, "one aspherical cycle: table, refine"),
  ("har_loop",   False, "full HAR loop, up to 4 cycles"),
  ("h_aniso",    False, "anisotropic hydrogens through one cycle"),
  ("h_free",     False, "hydrogens freed from AFIX through one cycle"),
  ("tsc_reuse",  False, "refine again from the table of the previous cycle"),
  ("no_refine",  False, "table only, run_refine off"),
  ("fcf_only",   False, "make_fcf_only: fcf from the table, no refinement"),
  ("alt_method", True,  "second method/basis of the same backend"),
  ("multipart",  True,  "PART disorder through -mtc (THPP)"),
]
DEFAULT_PATHS = "spherical,har1"


def _env_list(name, default):
  raw = os.environ.get(name, "")
  if not raw.strip():
    raw = default
  return [s.strip() for s in raw.split(",") if s.strip()]


def register(suite):
  full = os.environ.get("OLEX2_TEST_FULL", "") == "1"
  samples = _env_list("OLEX2_TEST_SAMPLES", DEFAULT_SAMPLES)
  wanted_backends = [s.lower() for s in _env_list("OLEX2_TEST_BACKENDS", "")]
  paths = _env_list("OLEX2_TEST_PATHS", DEFAULT_PATHS)
  if [p.lower() for p in paths] == ["all"]:
    paths = [p for p, slow, doc in PATHS]
  known_paths = [p for p, slow, doc in PATHS]
  for p in paths:
    if p not in known_paths:
      suite.run("matrix", "path %s" % p, _unknown_path, p, known_paths)
  paths = [p for p in paths if p in known_paths]

  suite.run("matrix", "NoSpherA2 present", base.t_present)
  suite.run("matrix", "experimental backends enabled", t_enable)
  suite.run("matrix", "backends offered", t_offered)

  for sample in samples:
    if sample not in SAMPLES:
      suite.run("matrix", "%s known" % sample, _unknown_sample, sample)
      continue
    if "spherical" in paths:
      suite.run("matrix", "%s spherical" % sample, t_spherical, suite, sample)
    for backend, slow_backend in BACKENDS:
      if wanted_backends and backend.lower() not in wanted_backends:
        continue
      for path in paths:
        if path in ("spherical", "multipart"):
          continue
        slow_path = dict((p, s) for p, s, d in PATHS)[path]
        suite.run("matrix", "%s %s %s" % (sample, backend, path),
                  t_path, suite, sample, backend, path,
                  slow_backend or slow_path, full)
    if "multipart" in paths and sample == "THPP":
      for backend in WAVEFUNCTION_BACKENDS:
        if wanted_backends and backend.lower() not in wanted_backends:
          continue
        suite.run("matrix", "%s %s multipart" % (sample, backend),
                  t_multipart, suite, sample, backend, full)


def _unknown_sample(sample):
  raise SkipTest("not a known sample (have %s)" % ", ".join(sorted(SAMPLES)))


def _unknown_path(path, known):
  raise AssertionError("no such path (have %s)" % ", ".join(known))


# ---------------------------------------------------------------------------
# environment
# ---------------------------------------------------------------------------

def t_enable():
  """Switch on what the GUI hides, for this session only."""
  done = []
  # SALTED and fragHAR are offered only in debug mode
  if not OV.IsDebugging():
    OV.SetParam('olex2.debug', True)
    done.append("olex2.debug")
  # OCC is added to the source list at start-up, from the user flag
  if OV.GetParam('user.NoSpherA2.show_OCC') != True:
    OV.SetParam('user.NoSpherA2.show_OCC', True)
    done.append("user.NoSpherA2.show_OCC")
  try:
    from NoSpherA2.NoSpherA2 import NoSpherA2_instance as nsp2
    nsp2.setup_occ_executables()
  except Exception as e:
    raise AssertionError("setup_occ_executables: %s" % e)
  model = os.environ.get("OLEX2_TEST_SALTED_MODEL", "").strip()
  if model:
    if not os.path.isdir(model):
      raise AssertionError("OLEX2_TEST_SALTED_MODEL is not a directory: %s" % model)
    listed = str(OV.GetParam('user.NoSpherA2.salted_models_list', '') or '')
    if model not in listed.replace(",", ";").split(";"):
      OV.SetParam('user.NoSpherA2.salted_models_list',
                  (listed + ";" + model) if listed.strip() else model)
    nsa2_set_param('selected_salted_model', model)
    done.append("SALTED model %s" % os.path.basename(model))
  ncpus = os.environ.get("OLEX2_TEST_NCPUS", "").strip()
  if not ncpus:
    try:
      import multiprocessing
      ncpus = str(min(8, multiprocessing.cpu_count()))
    except Exception:
      ncpus = "4"
  nsa2_set_param('ncpus', ncpus)
  nsa2_set_param('mem', os.environ.get("OLEX2_TEST_MEM", "4").strip() or "4")
  done.append("ncpus %s" % ncpus)
  return ", ".join(done) or "nothing to do"


def t_offered():
  offered = [s for s in _offered() if not s.startswith("Get ")]
  missing = [s[4:] for s in _offered() if s.startswith("Get ")]
  if not offered:
    raise AssertionError("NoSpherA2 offers no backend at all")
  return "%s%s" % (", ".join(offered),
                   (" (not found: %s)" % ", ".join(missing)) if missing else "")


# ---------------------------------------------------------------------------
# samples
# ---------------------------------------------------------------------------

def _sample_roots():
  # the directory the caller names wins over what the run directory ships:
  # the release data carries its own malbac (no solvent, an even electron
  # count) where sample_data's has the half-occupied toluene
  roots = []
  extra = os.environ.get("OLEX2_TEST_SAMPLE_DIR", "").strip()
  if extra:
    roots.append(extra)
  roots += [os.path.join(OV.BaseDir(), "sample_data"),
            os.path.join(OV.BaseDir(), "samples")]
  try:
    d = OV.GetParam('user.sample_dir', '')
    if d:
      roots.append(str(d))
  except Exception:
    pass
  return roots


def sample_copy(suite, name):
  """A fresh copy of a sample, looked up where olex2c can actually see one."""
  src = None
  for root in _sample_roots():
    cand = os.path.join(root, name)
    if os.path.isdir(cand):
      src = cand
      break
  if src is None:
    raise SkipTest("no sample %r under %s" % (name, "; ".join(_sample_roots())))
  dst = os.path.join(suite.scratch, name)
  if os.path.isdir(dst):
    shutil.rmtree(dst, ignore_errors=True)
  if os.path.isdir(dst):
    suite._copy_n += 1
    dst = "%s_%d" % (dst, suite._copy_n)
  # the movie folders of sucrose/timmy are not needed and not small
  shutil.copytree(src, dst, ignore=shutil.ignore_patterns("movie", "olex2", "*.tsc", "*.tscb"))
  return dst


def _model_file(folder):
  """res/ins first; a cif when that is all there is (epoxide)."""
  try:
    return model_in(folder)
  except SkipTest:
    for f in os.listdir(folder):
      if f.lower().endswith(".cif"):
        return os.path.join(folder, f)
    raise


def _load_model(folder, model):
  """Load and make the model refinable the way the GUI's refine macro does.

  The suite calls spy.refine.do_refine directly and so skips the macro's
  spy.AnalyseRefinementSource() step. For a res/ins that step is a no-op; for
  a sample that ships only a cif (epoxide) it writes the ins next to the cif
  and reloads it - without it Method.pre_refinement dies on olx.Ins('MORE')
  with 'INS file is expected' and no refinement runs at all.
  """
  load(model)
  macro("user '%s'" % folder.replace("\\", "/"))
  if model.lower().endswith(".cif"):
    import RunPrg
    if not RunPrg.AnalyseRefinementSource():
      raise AssertionError("%s could not be turned into a refinable ins"
                           % os.path.basename(model))
    if atom_count() == 0:
      raise AssertionError("no atoms after converting %s to ins"
                           % os.path.basename(model))


def _prepare(suite, sample):
  """Fresh copy, loaded, working directory set, refined spherically."""
  folder = sample_copy(suite, sample)
  if not has_hkl(folder):
    raise SkipTest("no hkl with the %s sample" % sample)
  model = _model_file(folder)
  _load_model(folder, model)
  # The GUI default carries the weighting scheme suggested by the previous
  # refinement into the next one (Method.py, update_weight). In one session
  # that is the previous *case*, so the spherical baseline of case N would
  # start from the aspherical weights of case N-1: cycle-1 wR2 drifted 6.6 ->
  # 7.0 across four sucrose copies of the same file. Every case refines with
  # the WGHT the sample ships instead.
  OV.SetParam('snum.refinement.update_weight', False)
  clear_r1()
  _refine_spherically()
  return folder, r1_of_last_refinement()


def _set_aspherical(sample, source):
  charge, mult, doc = SAMPLES[sample]
  nsa2_set_param('charge', str(charge))
  nsa2_set_param('multiplicity', str(mult))
  nsa2_set_param('source', source)
  # what the source combo of the GUI does on top of setting the source: the
  # chosen generator's own flags get defaults, every other generator's are
  # blanked. A tsc source (XCW's use button, say) blanks the ORCA SCF flags,
  # and ORCA's input writer needs them back before it can write a line.
  from utilities import reset_unused_generator_flags
  reset_unused_generator_flags(source)
  nsa2_set_param('full_HAR', False)
  nsa2_set_param('Max_HAR_Cycles', '1')
  nsa2_set_param('run_refine', True)
  nsa2_set_param('make_fcf_only', False)
  nsa2_set_param('h_aniso', False)
  nsa2_set_param('h_afix', False)
  nsa2_set_param('use_aspherical', True)


def _defaults_for(backend):
  """Method and basis as the GUI would leave them for a fresh structure."""
  if backend in ("ORCA", "OCC"):
    nsa2_set_param('method', 'B3LYP')
    nsa2_set_param('basis_name', 'def2-SVP')
  if backend == "ORCA":
    # the phil defaults, pinned so that the result does not depend on which
    # generator an earlier case left selected
    nsa2_set_param('ORCA_SCF_Conv', 'NoSpherA2SCF')
    nsa2_set_param('ORCA_SCF_Strategy', 'NormalConv')
    nsa2_set_param('ORCA_Solvation', 'Vacuum')
  elif backend == "xTB":
    nsa2_set_param('method', 'GFN2')


def _available(backend):
  """Resolve the backend name and what it needs, or skip saying why."""
  source = _source_name(backend)
  if backend == "SALTED":
    base._select_salted_model()
  elif backend == "OCC":
    share = os.path.join(OV.BaseDir(), "occ", "share", "basis")
    if not os.path.isdir(share):
      raise SkipTest("OCC data is not installed (%s)" % share)
  return source


# ---------------------------------------------------------------------------
# the paths
# ---------------------------------------------------------------------------

def t_spherical(suite, sample):
  folder, r1 = _prepare(suite, sample)
  if r1 > 0.20:
    raise AssertionError("spherical R1 %.4f on %s" % (r1, sample))
  return "R1 %.4f" % r1


def t_path(suite, sample, backend, path, slow, full):
  if slow and not full:
    raise SkipTest("quantum backend or slow path, set OLEX2_TEST_FULL=1")
  source = _available(backend)
  fn = globals()["path_" + path]
  return fn(suite, sample, backend, source)


def _aspherical_run(sample, backend, source, folder):
  """One deal_with_AAFF run with the current settings; the R1 it left."""
  clear_r1()
  said = _refine_capturing()
  cycles = _har_log_cycles(folder)
  if not cycles:
    raise AssertionError("%s: deal_with_AAFF left no HAR log%s"
                         % (backend, _reason(said)))
  tsc = _tsc_beside(folder)
  if tsc is None:
    raise AssertionError("%s produced no tsc or tscb%s" % (backend, _reason(said)))
  return cycles, tsc, said


def _check_moved(backend, r_sph, r_asp):
  if abs(r_asp - r_sph) < 1e-5:
    raise AssertionError("%s: R1 %.4f unchanged by the table - it did not reach "
                         "the refinement" % (backend, r_asp))
  if r_asp > r_sph + 0.01:
    raise AssertionError("%s: aspherical R1 %.4f is worse than spherical %.4f"
                         % (backend, r_asp, r_sph))


def path_har1(suite, sample, backend, source):
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_asp)
  return "%s, %d cycle, R1 %.4f -> %.4f" % (os.path.basename(tsc), cycles, r_sph, r_asp)


def path_har_loop(suite, sample, backend, source):
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  nsa2_set_param('full_HAR', True)
  nsa2_set_param('Max_HAR_Cycles', '4')
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_asp)
  log = _har_log_text(folder)
  if cycles > 4:
    raise AssertionError("%d cycles recorded for Max_HAR_Cycles 4" % cycles)
  unconverged = "UNCONVERGED" in log
  if cycles < 4 and unconverged:
    raise AssertionError("stopped after %d cycle(s) without converging" % cycles)
  return "%d cycle(s)%s, R1 %.4f -> %.4f" % (
    cycles, " (unconverged at the cap)" if unconverged else " converged", r_sph, r_asp)


def path_h_aniso(suite, sample, backend, source):
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  nsa2_set_param('h_aniso', True)
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_asp)
  n_h, n_aniso = _hydrogens()
  if n_h and n_aniso < n_h:
    raise AssertionError("%d of %d hydrogens are still isotropic after h_aniso"
                         % (n_h - n_aniso, n_h))
  return "%d/%d H anisotropic, R1 %.4f -> %.4f" % (n_aniso, n_h, r_sph, r_asp)


def path_h_free(suite, sample, backend, source):
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  nsa2_set_param('h_afix', True)
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_asp)
  riding = _riding_hydrogens()
  if riding:
    raise AssertionError("%d hydrogens still carry an AFIX after h_afix" % riding)
  return "no AFIX H left, R1 %.4f -> %.4f" % (r_sph, r_asp)


def path_tsc_reuse(suite, sample, backend, source):
  """The table of a finished run, chosen as the source, refines again."""
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_first = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_first)
  # now the way a user picks a table from the drop-down: the source *is* the file
  table = os.path.basename(tsc)
  nsa2_set_param('source', table)
  nsa2_set_param('file', table)
  clear_r1()
  said = _refine_capturing()
  r_second = r1_of_last_refinement()
  if _har_log_cycles(folder) < 1:
    raise AssertionError("no HAR log from the reused table%s" % _reason(said))
  if abs(r_second - r_first) > 0.005:
    raise AssertionError("R1 from the reused table %.4f differs from the run "
                         "that wrote it %.4f" % (r_second, r_first))
  return "%s reused, R1 %.4f -> %.4f -> %.4f" % (table, r_sph, r_first, r_second)


def path_no_refine(suite, sample, backend, source):
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  nsa2_set_param('run_refine', False)
  clear_r1()
  said = _refine_capturing()
  tsc = _tsc_beside(folder)
  if tsc is None:
    raise AssertionError("%s produced no table%s" % (backend, _reason(said)))
  r_after = _r1_raw()
  if r_after is not None and abs(r_after - r_sph) > 1e-6:
    raise AssertionError("run_refine off, but R1 moved %.4f -> %.4f" % (r_sph, r_after))
  return "%s written, no refinement" % os.path.basename(tsc)


def path_fcf_only(suite, sample, backend, source):
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  _defaults_for(backend)
  # a table first, as the fcf route reads the table it finds
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  name = os.path.splitext(os.path.basename(_model_file(folder)))[0]
  for f in os.listdir(folder):
    if f.lower().endswith(".fcf"):
      os.remove(os.path.join(folder, f))
  nsa2_set_param('make_fcf_only', True)
  try:
    said = _refine_capturing()
  finally:
    nsa2_set_param('make_fcf_only', False)
  fcfs = [f for f in os.listdir(folder) if f.lower().endswith(".fcf")]
  if not fcfs:
    raise AssertionError("make_fcf_only wrote no fcf%s" % _reason(said))
  n = 0
  for line in open(os.path.join(folder, fcfs[0]), "r", errors="ignore"):
    parts = line.split()
    if len(parts) >= 5 and parts[0].lstrip("-").isdigit():
      n += 1
  if n == 0:
    raise AssertionError("%s carries no reflections" % fcfs[0])
  return "%s with %d reflections (R1 %.4f)" % (fcfs[0], n, r_asp)


def path_alt_method(suite, sample, backend, source):
  if backend not in ALT_METHOD:
    raise SkipTest("%s has no second method to run" % backend)
  method, basis = ALT_METHOD[backend]
  folder, r_sph = _prepare(suite, sample)
  _set_aspherical(sample, source)
  nsa2_set_param('method', method)
  if basis:
    nsa2_set_param('basis_name', basis)
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  r_asp = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_asp)
  return "%s%s, R1 %.4f -> %.4f" % (method, ("/" + basis) if basis else "", r_sph, r_asp)


def t_multipart(suite, sample, backend, full):
  if not full:
    raise SkipTest("needs a wavefunction backend, set OLEX2_TEST_FULL=1")
  source = _available(backend)
  folder = sample_copy(suite, sample)
  if not has_hkl(folder):
    raise SkipTest("no hkl with %s" % sample)
  model = _model_file(folder)
  parts = set()
  for line in open(model, "r", errors="ignore"):
    if line[:5].strip().upper() == "PART":
      parts.add(line.split()[1])
  if len(parts) < 2:
    raise SkipTest("%s no longer carries more than one part" % sample)
  _load_model(folder, model)
  clear_r1()
  _refine_spherically()
  r_sph = r1_of_last_refinement()
  _set_aspherical(sample, source)
  _defaults_for(backend)
  cycles, tsc, said = _aspherical_run(sample, backend, source, folder)
  leftovers = _per_part_tables(folder)
  if leftovers:
    raise AssertionError("per-part tables were left to be merged: %s"
                         % ", ".join(sorted(os.path.basename(p) for p in leftovers)))
  wfns = _wavefunctions(folder)
  if len(wfns) < 2:
    raise AssertionError("a %d-part structure produced %d wavefunction(s)"
                         % (len(parts), len(wfns)))
  tables = _all_tables(folder)
  if len(tables) != 1:
    raise AssertionError("expected one table, found %d: %s"
                         % (len(tables), ", ".join(sorted(os.path.basename(t) for t in tables))))
  r_asp = r1_of_last_refinement()
  _check_moved(backend, r_sph, r_asp)
  return "parts %s -> one %s, R1 %.4f -> %.4f" % (
    ",".join(sorted(parts)), os.path.basename(tsc), r_sph, r_asp)


# ---------------------------------------------------------------------------
# model queries
# ---------------------------------------------------------------------------

def _har_log_text(folder):
  text = []
  for root, dirs, files in os.walk(folder):
    for f in files:
      if f.endswith(".NoSpherA2"):
        text.append(open(os.path.join(root, f), "r", errors="ignore").read())
  return "\n".join(text)


def _live_atoms():
  n = int(olx.xf.au.GetAtomCount())
  for i in range(n):
    if str(olx.xf.au.IsAtomDeleted(i)).lower() == "true":
      continue
    if str(olx.xf.au.IsPeak(i)).lower() == "true":
      continue
    yield i


def _hydrogens():
  """(hydrogens, of which anisotropic) - GetAtomU is one value or six."""
  n_h = n_aniso = 0
  for i in _live_atoms():
    if str(olx.xf.au.GetAtomType(i)).strip() not in ("H", "D"):
      continue
    n_h += 1
    u = str(olx.xf.au.GetAtomU(i))
    if u.count(",") >= 5:
      n_aniso += 1
  return n_h, n_aniso


def _riding_hydrogens():
  n = 0
  for i in _live_atoms():
    if str(olx.xf.au.GetAtomType(i)).strip() not in ("H", "D"):
      continue
    try:
      afix = int(str(olx.xf.au.GetAtomAfix(i)))
    except ValueError:
      afix = 0
    if afix:
      n += 1
  return n
