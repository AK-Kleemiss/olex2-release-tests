"""Pictures of the 3D view: the model, its maps and the NoSpherA2 surfaces.

The workflow case goes the whole way first - FLINT, hydrogens, a spherical
refinement, a wavefunction, one HAR cycle, the properties as cubes - so that the
maps it then photographs are the ones a run computes rather than the ones a
sample ships.

`matr` sets the camera - from the cell (matr a/b/c) or, in its nine-number
form, from the Cartesian axes, which does not depend on the cell at all - so
every run frames the same thing and two runs of the suite can be held side by
side. `pict` then renders the GL canvas alone, without the panels and the
console, at a fixed width into OLEX2_TEST_GUI_DIR, where the driver's
index.html lists it next to the GUI screenshots.

The golden holds what a script can decide without a reference picture: that
every render arrived and is not an empty canvas, how many of the cameras gave
a different picture (one that never turned makes two of them byte-identical),
that the surface or the map changed the picture at all, and the numbers behind
the scene - the faces of a surface. Whether the picture looks *right* - a hole
in a surface, an inverted normal, a colour ramp the wrong way round - is for
eyes, and that is what the pictures are for.

Without a GUI (olex2c) every case skips: pict needs a GL context.
"""
from __future__ import absolute_import, division, print_function

import os
import shutil
import sys

import olx
from olexFunctions import OV

from pipeline_tests import SkipTest
from gui_tests import gui_dir, log_size, log_errors, _assert_clean, snapshot

# pict reads a second number above 100 as the target width, so the size of the
# Olex2 window does not decide the picture
PICT_WIDTH = "900"
# at that width a canvas holding nothing but the background colour compresses
# to about 2 kB; anything drawn on it is far more
MIN_PNG = 8000

# down each cell axis, and one camera that is the same on every structure
VIEWS = [("a", ("a",)), ("b", ("b",)), ("c", ("c",)),
         ("cart", ("1", "0", "0", "0", "1", "0", "0", "0", "1"))]


def _load(suite, sample):
  """A GUI and a known structure on screen with an empty scene: what is
  loaded decides the picture, so no case here draws whatever the previous
  one happened to leave behind."""
  if not OV.HasGUI():
    raise SkipTest("no GUI (console build)")
  olx.app.AutoFlushLog("true")
  from group_nsa2_matrix import sample_copy, _model_file, _load_model
  folder = sample_copy(suite, sample)
  _load_model(folder, _model_file(folder))
  os.chdir(folder)
  import cubes_maps
  cubes_maps.disable_map()
  return folder


def render(state, view=None):
  """Orient with matr, render the canvas into the gui dir, hand back the file.

  The bytes are the return value on purpose: wxWidgets writes the png
  deterministically, so two states that produced the same bytes drew the same
  thing - which is how a camera that never turned, or a surface that never
  reached the scene, is caught without a reference picture."""
  path = os.path.join(gui_dir(), state + ".png")
  if os.path.isfile(path):
    os.remove(path)
  if view:
    olx.Matr(*view)
  olx.Refresh()
  if sys.platform == "win32":
    olx.Pict(path, PICT_WIDTH)  # the GL canvas alone; the macro is Windows-only
  if not os.path.isfile(path):
    snapshot(state)  # elsewhere: the window shot, panels and all
  if not os.path.isfile(path):
    raise AssertionError("no picture written for %s" % state)
  size = os.path.getsize(path)
  if size < MIN_PNG:
    raise AssertionError("%s is %d bytes - the canvas came out empty" % (state, size))
  with open(path, "rb") as f:
    return f.read()


def _cameras(prefix, views=VIEWS):
  return [render("%s_%s" % (prefix, name), args) for name, args in views]


def _seen(views, pics, *extra):
  pairs = list(extra) + [("views", ",".join(n for n, _ in views)),
                         ("distinct", len(set(pics))), ("errors", 0)]
  return " ".join("%s=%s" % (k, v) for k, v in pairs)


# ---------------------------------------------------------------------------
# cases
# ---------------------------------------------------------------------------

def c_render_views(suite):
  """The model itself, down each cell axis and from the Cartesian camera."""
  _load(suite, "epoxide")
  at = log_size()
  pics = _cameras("view")
  _assert_clean(log_errors(at), "model views")
  return _seen(VIEWS, pics)


def c_render_hirshfeld(suite):
  """The Hirshfeld surface of the molecule against its crystal environment,
  coloured by d_norm. Promolecular, so no wavefunction is needed."""
  _load(suite, "epoxide")
  from group_release import _params
  import cubes_maps
  at = log_size()
  bare = render("hirshfeld_bare", VIEWS[0][1])  # the model alone, same camera
  with _params(**{"map.resolution": "0.4"}):
    cubes_maps.hirshfeld_surface("Hirshfeld d_norm")
    values = cubes_maps._hs_values("Hirshfeld d_norm")
    if values is None or not len(values):
      raise AssertionError("the run left no Hirshfeld_surface.dat with a d_norm column")
    pics = _cameras("hirshfeld", VIEWS[:3])
  if pics[0] == bare:
    raise AssertionError("the picture is what it was without the surface - "
                         "it never reached the scene")
  _assert_clean(log_errors(at), "Hirshfeld surface")
  return _seen(VIEWS[:3], pics, ("faces", len(values)))


def c_render_esp_surface(suite):
  """The electron density isosurface coloured by the electrostatic potential:
  the one surface that needs a density source - a wavefunction beside the
  structure, or the SALTED model the run was given."""
  _load(suite, "epoxide")
  from group_release import _params
  import cubes_maps
  if not cubes_maps.has_density_source():
    raise SkipTest("no wavefunction and no SALTED model selected")
  at = log_size()
  bare = render("esp_surface_bare", VIEWS[0][1])
  with _params(**{"map.type": "ESP surface", "map.resolution": "0.4"}):
    cubes_maps.change_map()
    obj = cubes_maps._surface_obj("ESP surface")
    if not os.path.isfile(obj):
      raise AssertionError("no %s written" % os.path.basename(obj))
    faces = sum(1 for line in open(obj, errors="ignore") if line.startswith("f "))
    pics = _cameras("esp_surface", VIEWS[:3])
  if pics[0] == bare:
    raise AssertionError("the picture is what it was without the surface - "
                         "it never reached the scene")
  _assert_clean(log_errors(at), "ESP surface")
  return _seen(VIEWS[:3], pics, ("faces", faces))


def c_render_map(suite):
  """The residual density grid after a spherical refinement, the way the map
  controls of the GUI put it on screen."""
  if not OV.HasGUI():
    raise SkipTest("no GUI (console build)")
  olx.app.AutoFlushLog("true")
  from group_nsa2_matrix import _prepare
  import cubes_maps
  _prepare(suite, "epoxide")
  cubes_maps.disable_map()
  at = log_size()
  bare = render("map_bare", VIEWS[0][1])
  cubes_maps.show_fft_map(0.3, "diff")
  pics = _cameras("map_diff", VIEWS[:3])
  if pics[0] == bare:
    raise AssertionError("the picture is what it was without the map - "
                         "nothing was put on screen")
  _assert_clean(log_errors(at), "residual map")
  return _seen(VIEWS[:3], pics, ("map", "diff"))


# ---------------------------------------------------------------------------
# the whole workflow, then a picture of every map it computed
# ---------------------------------------------------------------------------

# the cubes a Property_ flag asks NoSpherA2 for; -rdg writes the signed density
# beside the gradient, and it is the pair of them that makes NCI appear
PROPERTY_CUBES = [("Property_Lap", ["lap"]), ("Property_Eli", ["eli"]),
                  ("Property_Elf", ["elf"]), ("Property_RDG", ["rdg", "signed_rho"]),
                  ("Property_ESP", ["esp"]), ("Property_DEF", ["def"]),
                  ("Property_MO", ["MO_0"])]
# SALTED predicts a density, not orbitals, and calculate_cubes turns these back
# off with a line in the log - so a SALTED workflow can only reach the other three
NEEDS_WFN = ("Property_Elf", "Property_RDG", "Property_DEF", "Property_MO")
WORKFLOW_CUBE_TIMEOUT = 900.0   # seconds for the seven cubes of a properties run
# the surfaces have their own cases above, and the intermolecular NCI is a run
# over the whole displayed cluster - more than a picture is worth here
SKIP_TYPES = ("Intermolecular NCI",)


def _slug(text):
  return "_".join("".join(c if c.isalnum() else " " for c in text).split()).lower()


def _wfn_beside(folder):
  """cubes_maps looks for the wavefunction beside the structure, the generator
  leaves it in its job folder: put a copy where it is looked for."""
  import cubes_maps
  if cubes_maps._wfn_file():
    return True
  name = OV.ModelSrc()
  for ext in (".gbw", ".wfx", ".ffn", ".wfn", ".molden", ".xtb"):
    src = os.path.join(folder, "olex2", "Wfn_job", name + ext)
    if os.path.isfile(src):
      shutil.copy(src, os.path.join(folder, name + ext))
      return True
  return False


def c_render_workflow(suite, sample, backend):
  """From the hkl to the pictures: FLINT, hydrogens, a spherical
  refinement, the wavefunction of `backend`, one HAR cycle, the properties as
  cubes - and then one photograph of every map type the run put on the menu.

  None of these maps is shipped with the sample: the ELI-D, the Laplacian and
  their neighbours come out of the properties run, the residual and its four
  siblings are transformed from the refinement each time they are shown, and
  Rho + ESP is computed on the spot when it is picked. The golden pins which
  of them the workflow reached and that each one changed the picture; whether
  the ELI-D basins sit where they belong is for eyes, and that is what the
  pictures are for.
  """
  if not OV.HasGUI():
    raise SkipTest("no GUI (console build)")
  olx.app.AutoFlushLog("true")
  import group_flint
  import group_nsa2_matrix as matrix
  import cubes_maps
  from group_release import _params, _settle, _check_not_worse
  from group_nsa2_matrix import _set_aspherical, _defaults_for, _aspherical_run
  from group_nosphera2 import _refine_spherically
  from pipeline_tests import macro, atom_count, clear_r1, r1_of_last_refinement

  prefix = "workflow_" + backend.lower()
  source = matrix._available(backend)
  at = log_size()

  # solve - the case of its own that owns the numbers; here it is a step
  group_flint.t_flint(suite, sample)
  folder = OV.FilePath()
  os.chdir(folder)
  solved = atom_count()
  macro("HAdd")   # a solution has no peaks for the hydrogens; the HAR needs them
  hydrogens = atom_count() - solved

  OV.SetParam('snum.refinement.update_weight', False)
  clear_r1()
  _refine_spherically()
  r_sph = r1_of_last_refinement()
  _set_aspherical(sample, source)
  _defaults_for(backend)
  _aspherical_run(sample, backend, source, folder)
  # a solved model is not the deposited one, so R1 is only held to not getting
  # worse - that the table reached the refinement is what _aspherical_run checks
  _check_not_worse(backend, r_sph, r1_of_last_refinement())

  name = OV.ModelSrc()
  wfn = _wfn_beside(folder)
  flags = dict((k, bool(wfn) or k not in NEEDS_WFN) for k, _ in PROPERTY_CUBES)
  wanted = [c for k, cs in PROPERTY_CUBES if flags[k] for c in cs]
  cubes = [cubes_maps._p("%s_%s.cube" % (name, w)) for w in wanted]
  for c in cubes:   # a cube of an earlier case must not pass as this run's
    for f in (c, c + "b"):
      if os.path.isfile(f):
        os.remove(f)
  flags.update({"Property_ATOM": False, "Property_all_MOs": False,
                "Property_ESP_surface": False,
                "map.radius": "1.0", "map.resolution": "0.5"})
  with _params(**flags):
    cubes_maps.calculate_cubes()
  _settle(cubes, timeout=WORKFLOW_CUBE_TIMEOUT)

  skip = set(cubes_maps.SURFACE_OBJS) | set(SKIP_TYPES)
  offered = [t for t in str(cubes_maps.get_map_types()).split(";") if t and t not in skip]
  if not offered:
    raise AssertionError("the properties run left no map type on the menu")
  cubes_maps.disable_map()
  bare = render(prefix + "_bare", VIEWS[0][1])
  pics, shown = [], []
  for entry in offered:
    # the combo is built as "Residual<-diff": the name is shown, the value is set
    label, _, value = entry.partition("<-")
    with _params(**{"map.type": value or label, "map.resolution": "0.5"}):
      cubes_maps.change_map()
      pic = render("%s_%s" % (prefix, _slug(label)), VIEWS[0][1])
    if pic == bare:
      raise AssertionError("%s: the picture is what it was with no map - "
                           "it never reached the scene" % label)
    pics.append(pic)
    shown.append(label)
  _assert_clean(log_errors(at), "%s workflow" % backend)
  return " ".join(["source=%s" % backend, "solved=%d" % solved,
                   "hydrogens=%d" % hydrogens, "cubes=%s" % ",".join(wanted),
                   "maps=%s" % ",".join(shown),
                   "distinct=%d" % len(set(pics)), "errors=0"])


# name, tier, function, keyword arguments - the order is the golden's order
CASES = [
  ("render_views",           "quick", c_render_views, {}),
  ("render_hirshfeld",       "quick", c_render_hirshfeld, {}),
  ("render_esp_surface",     "full",  c_render_esp_surface, {}),
  ("render_map",             "full",  c_render_map, {}),
  ("render_workflow_salted", "full",  c_render_workflow,
   {"sample": "epoxide", "backend": "SALTED"}),
  ("render_workflow_orca",   "full",  c_render_workflow,
   {"sample": "epoxide", "backend": "ORCA"}),
]
