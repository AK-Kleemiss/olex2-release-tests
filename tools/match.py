import sys
from iotbx import shelx
from cctbx.eltbx import tiny_pse
def read(fn):
  if fn.lower().endswith(".cif"):
    from iotbx import cif
    return list(cif.reader(file_path=fn).build_crystal_structures().values())[0]
  return shelx._cctbx_xray_structure_from(filename=fn)
def Z(t):
  try: return tiny_pse.table(t.strip().capitalize()).atomic_number()
  except Exception: return 0
ref, got = read(sys.argv[1]), read(sys.argv[2])
uc, sg = ref.unit_cell(), ref.space_group()
refs = [(s.label, s.scattering_type.strip(), s.site) for s in ref.scatterers() if s.scattering_type.strip() not in ("H", "D")]
sol = [(s.label, s.scattering_type.strip(), s.site, s.u_iso) for s in got.scatterers() if s.scattering_type.strip() != "Q"]
def wrap(v, w): return [a - round(a - b) for a, b in zip(v, w)]
def nearest(site, shift, inv):
  best = (9e9, None)
  for lab, el, rs in refs:
    for op in sg:
      x = op*rs
      if inv: x = [-a for a in x]
      x = [a + b for a, b in zip(x, shift)]
      d = uc.distance(wrap(x, site), site)
      if d < best[0]: best = (d, (lab, el))
  return best
def score(shift, inv, subset=None):
  return sum(1 for l, t, s, u in (subset or sol) if nearest(s, shift, inv)[0] < 0.5)
probe = sol[::max(1, len(sol)//8)][:8]
heavy = sorted(sol, key=lambda a: -Z(a[1]))[:3]
best = (-1, (0, 0, 0), False)
seen = set()
for inv in (False, True):
  for l, t, s, u in heavy:
    for lab, el, rs in refs:
      for op in sg:
        x = op*rs
        if inv: x = [-a for a in x]
        shift = tuple(round(a - b, 3) for a, b in zip(s, x))
        if (shift, inv) in seen: continue
        seen.add((shift, inv))
        if score(shift, inv, probe) < len(probe)//2: continue
        sc = score(shift, inv)
        if sc > best[0]: best = (sc, shift, inv)
sc, shift, inv = best
print("origin shift %s inverted=%s matched %d/%d" % (shift, inv, sc, len(sol)))
for l, t, s, u in sol:
  d, (lab, el) = nearest(s, shift, inv)
  flag = "" if el == t and d < 0.5 else ("  <-- %s" % el if d < 0.5 else "  <-- NOISE")
  print("%-5s %-3s U=%7.4f  nearest %-5s %-3s %.2f A%s" % (l, t, u, lab, el, d, flag))
