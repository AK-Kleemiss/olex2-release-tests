"""FLINT on the sample structures: the multi-trial charge-flipping
pipeline with its space-group shortlist and element assignment.

Each case loads the deposited model as the reference, runs
`olex2.solve` / `FLINT` through spy.RunSolutionPrg() exactly as the GUI
does, lets the deferred tidy-up (compaq, four cycles, ADP prune, re-typing)
run, and then compares what came back with what was deposited:

  sg_rank      where the deposited space group sits in the shortlist (0 = absent)
  cc           best correlation over the trials (seeds are fixed, so stable)
  peaks        peaks the solution posted, and how many of them the geometry
               classifier moved away from the density call
  allowed      the formula restriction the assignment ran under
  pruned       peaks the ADP prune removed
  retyped      labels the post-cleanup re-typing changed (none: the ADP
               prune never settled within its rounds, so no re-typing ran)
  types/typed  element histogram of the final model, and the fraction of the
               deposited non-H atoms it accounts for

FLINT is opt-in (user.solution.flint) and needs a cctbx with
smtbx.ab_initio; without that the group skips. The geometry classifier needs
etc/geometry_aid_model.npz beside NoSpherA2.exe, and a run that falls back to
density only is a failure, not a variant.

OLEX2_TEST_FLINT_SAMPLES  comma-separated sample names, default SAMPLES
                              + INORGANIC; "inorganic" = INORGANIC alone

The inorganic subsection (INORGANIC, COD entries, 21 Sep 2026) is the
carbon-free path: no C-C pair sets the density scale, so the formula ranks
do, the refined U median sits far under the organic 0.03-0.06 and heavy
atoms typed light collapse it further, and the re-typing arbitrates between
neighbouring Z (Na/Al/Si, Ca/Te) without a geometry class. Six cases are
the loss mechanisms seen on 3000 COD entries, two are all-right guards:

  cod_1560875  Gd O14 P5             prune cascade (median U -> 0), Gd<->P
  cod_2241658  Ag0.6 Fe Mo2 Na0.4 O8  mixed site, O->Na->Fe drift
  cod_2108240  plagioclase Al/Si/Na/Ca  adjacent Z, I-1 setting
  cod_2013004  K2 Mn2 O21 P6 Sr3      44 peaks pruned to 14, P21 at rank 3
  cod_2104335  Ca O3 Te              Ca<->Te swaps, P43 at rank 2
  cod_2208447  As4 Cs4 Se8           negative median U, Se->Cs
  cod_2229150  H O5 Pr S             all right (guard)
  cod_2108989  B4 Bi0.07 Fe3 O12 Sm0.93  all right but B (guard, R32)
"""
from __future__ import absolute_import, division, print_function

import os
import re

import olex
import olx
from olexFunctions import OV

from pipeline_tests import (macro, SkipTest, atom_count, space_group, has_hkl)
from group_nsa2_matrix import sample_copy, _model_file, _load_model
from group_nosphera2 import _refine_capturing

GROUP = "flint"
SAMPLES = ("sucrose", "epoxide", "water", "malbac")
INORGANIC = ("cod_1560875", "cod_2241658", "cod_2108240", "cod_2013004",
             "cod_2104335", "cod_2208447", "cod_2229150", "cod_2108989")
NPZ = os.path.join("etc", "geometry_aid_model.npz")


def register(suite):
  wanted = os.environ.get("OLEX2_TEST_FLINT_SAMPLES", "").strip()
  samples = [s.strip() for s in wanted.split(",") if s.strip()]
  samples = INORGANIC if samples == ["inorganic"] else samples or SAMPLES + INORGANIC
  for s in samples:
    suite.run(GROUP, "flint %s" % s, t_flint, suite, s)


def _method():
  """The FLINT method, registered for this session."""
  try:
    import smtbx.ab_initio  # noqa: F401
  except ImportError:
    raise SkipTest("this cctbx has no smtbx.ab_initio")
  import ExternalPrgParameters as EPP
  OV.SetParam('user.solution.flint', True)
  EPP.SPD, EPP.RPD = EPP.defineExternalPrograms()
  prg = EPP.SPD.programs.get('olex2.solve')
  method = prg.methods.get('FLINT') if prg else None
  if method is None:
    raise AssertionError("FLINT is not registered with olex2.solve")
  return method


def _types():
  """element -> atoms of the model, H and Q peaks left out."""
  out = {}
  for i in range(atom_count()):
    t = str(olx.xf.au.GetAtomType(i))
    if t in ("H", "D", "Q"):
      continue
    out[t] = out.get(t, 0) + 1
  return out


def _histogram(types):
  return ",".join("%s:%d" % (k, types[k]) for k in sorted(types))


def _grab(pattern, text, default=None, cast=int):
  m = re.search(pattern, text)
  return cast(m.group(1)) if m else default


def _expand_hkl_to_p1():
  # the CIF-derived hkl holds one asymmetric unit; a P1 solve needs the sphere
  from iotbx.shelx import hklf
  from cctbx import crystal
  src = OV.HKLSrc()
  cs = crystal.symmetry(
    unit_cell=[float(x) for x in olex.f("xf.au.GetCell()").split(',')],
    symbol="hall: " + olex.f("sg(%HS)"))
  ma = hklf.reader(file_name=src).as_miller_arrays(
    crystal_symmetry=cs)[0]
  with open(src, "w") as f:
    ma.expand_to_p1().export_as_shelx_hklf(f)


def t_flint(suite, sample):
  method = _method()
  if not os.path.isfile(os.path.join(OV.BaseDir(), NPZ)):
    raise AssertionError("%s is not in the run directory" % NPZ)
  folder = sample_copy(suite, sample)
  if not has_hkl(folder):
    raise SkipTest("no hkl with the %s sample" % sample)
  _load_model(folder, _model_file(folder))
  ref_types, ref_sg = _types(), space_group()
  ref_no = int(olex.f("sg(%#)"))
  if ref_no < 1:
    # a setting outside Olex2's table (I -1) still has a group type
    from cctbx import sgtbx
    ref_no = sgtbx.space_group_info(
      symbol="hall: " + olex.f("sg(%HS)")).type().number()
  n_ref = sum(ref_types.values())
  # OLEX2_TEST_FLINT_SG=P1 solves in a lower group to measure whether the
  # deposited one is recovered into the shortlist
  lower = os.environ.get("OLEX2_TEST_FLINT_SG", "").strip()
  if lower:
    _expand_hkl_to_p1()
    macro("ChangeSG %s" % lower)

  OV.SetParam('snum.solution.program', 'olex2.solve')
  OV.SetParam('snum.solution.method', 'FLINT')
  OV.SetParam('snum.solution.retype_after_tidy', True)
  # The GUI asks "solve this again?" once a model is loaded; that would block
  # a headless run in the message box.
  ask = OV.GetParam('user.alert_solve_anyway')
  OV.SetParam('user.alert_solve_anyway', 'N')
  try:
    text = _refine_capturing("spy.RunSolutionPrg()")
  finally:
    OV.SetParam('user.alert_solve_anyway', ask)

  if "No solution found" in text or atom_count() == 0:
    raise AssertionError("FLINT found no solution for %s" % sample)
  trials = _grab(r"Best of (\d+) trial", text)
  if not trials or (trials < 2 and "good_enough" not in text):
    raise AssertionError("expected several trials, log says %r" % trials)
  cc = _grab(r"Best of \d+ trial\(s\): correlation ([0-9.]+)", text, cast=float)
  if "using density only" in text or "Geometry step unavailable" in text:
    raise AssertionError("element assignment fell back to density only")
  peaks = _grab(r"Element assignment: (\d+) peaks", text)
  if peaks is None:
    raise AssertionError("no element assignment ran")
  changed = _grab(r"Element assignment: \d+ peaks, (\d+) where geometry", text, 0)
  allowed = _grab(r"Element assignment restricted to: ([^\n]+)", text, "any",
                  cast=lambda s: ",".join(s.replace(",", " ").split()))
  pruned = _grab(r"Pruned (\d+) peak", text, 0)
  if "Re-typed" in text:
    retyped = _grab(r"(\d+) label\(s\) changed", text, 0)
  elif "Re-typing skipped" in text:
    retyped = "skipped"
  else:
    retyped = "none"

  solver = method.cctbx_solver
  sugg = getattr(solver, 'solution_suggestions', None)
  entries = list(getattr(sugg, 'suggestions', []) or [])
  if not entries:
    raise AssertionError("no space-group suggestions")
  numbers = [e.space_group_info.type().number() for e in entries]
  sg_rank = numbers.index(ref_no) + 1 if ref_no in numbers else 0
  if not sg_rank:
    raise AssertionError("deposited %s (No. %d) is not among the suggestions %s"
                         % (ref_sg, ref_no, numbers))

  # The chooser table is built only with a GUI, so build it here: every
  # candidate placed, R1 column present (21 Sep 2026: R1 replaced the peak count).
  chooser = "-"
  if len(entries) >= 2:
    written = solver.writeSuggestions(solver.solution_f_obs, sugg)
    html = solver.suggestionsTableHtml(written, sugg)
    if not written or "<b>R1</b>" not in html:
      raise AssertionError("chooser table: %d candidates placed, R1 %s"
                           % (len(written), "<b>R1</b>" in html))
    # a candidate the P1 solution cannot be symmetrised into is left out
    chooser = ",".join("%.3f" % r1 for _, _, r1 in written)

  addsym = _grab(r"Possible missed symmetry: the refined model has ([^;]+), "
                 r"solved in", text, "-", cast=lambda s: s.replace(" ", ""))
  got = _types()
  typed = sum(min(got.get(k, 0), ref_types[k]) for k in ref_types)/float(n_ref)
  if len(ref_types) > 1 and list(got) == ["C"]:
    raise AssertionError("every peak stayed carbon; the formula is %s"
                         % _histogram(ref_types))
  if typed < 0.5:
    raise AssertionError("final model %s accounts for %.2f of the deposited %s"
                         % (_histogram(got), typed, _histogram(ref_types)))
  return " ".join("%s=%s" % kv for kv in [
    ("sg", ref_sg), ("sg_rank", sg_rank), ("n_suggest", len(entries)),
    ("chooser_r1", chooser),
    ("addsym", addsym),
    ("trials", trials), ("cc", "%.3f" % cc), ("peaks", peaks),
    ("geometry_changed", changed), ("allowed", allowed), ("pruned", pruned),
    ("retyped", retyped), ("atoms", sum(got.values())),
    ("types", _histogram(got)), ("typed", "%.2f" % typed)])
