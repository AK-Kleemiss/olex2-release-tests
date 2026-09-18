# stage_cod.py <split-file> <n> <out-ids>: local COD mirror entries -> D:/devel/olex2-samples/cod_<id>
import random, re, sys, os, shutil, subprocess
random.seed(int(sys.argv[4]) if len(sys.argv) > 4 else 17)
COD = 'C:/Users/florian/cctbx-work/COD'
ids = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
skip = set(l.split()[0] for f in sys.argv[5:] for l in open(f) if l.strip())
random.shuffle(ids)
want, kept = int(sys.argv[2]), []
def path(kind, cid): return '%s/%s/%s/%s/%s/%s.%s' % (COD, kind, cid[0], cid[1:3], cid[3:5], cid, kind)
for cid in ids:
  if len(kept) >= want: break
  if cid in skip or not os.path.isfile(path('hkl', cid)) or not os.path.isfile(path('cif', cid)): continue
  cif = open(path('cif', cid), encoding='utf-8', errors='replace').read()
  if cif.count('\ndata_') != 1: continue
  if not re.search(r'_diffrn_radiation_type\s+[\'"]?(Mo|Cu)', cif): continue
  m = re.search(r'_refine_ls_R_factor_gt\s+([0-9.]+)', cif)
  if not m or float(m.group(1)) > 0.08: continue
  lines = cif.splitlines()
  try: i = next(k for k, l in enumerate(lines) if l.startswith('_atom_site_label'))
  except StopIteration: continue
  while lines[i].startswith('_atom_site_'): i += 1
  rows = []
  while i < len(lines) and lines[i].strip() and not lines[i].startswith(('_', 'loop_', '#')):
    rows.append(lines[i].split()); i += 1
  nonh = [r for r in rows if len(r) > 1 and r[1].strip("'").capitalize() not in ('H', 'D')]
  if not 8 <= len(nonh) <= 70: continue
  hkl = open(path('hkl', cid), encoding='utf-8', errors='replace').read()
  if '_refln_F_squared_meas' not in hkl or '_refln_F_squared_sigma' not in hkl: continue
  d = 'D:/devel/olex2-samples/cod_%s' % cid
  os.makedirs(d, exist_ok=True)
  open('%s/%s.cif' % (d, cid), 'w', encoding='utf-8').write(cif)
  open('%s/%s.hkl' % (d, cid), 'w', encoding='utf-8').write(hkl)
  r = subprocess.run([sys.executable, 'cod2hkl.py', '%s/%s.hkl' % (d, cid)], capture_output=True, text=True)
  if r.returncode: shutil.rmtree(d); print('drop', cid, r.stderr.strip()[-80:]); continue
  kept.append(cid); print(cid, len(nonh), m.group(1), flush=True)
open(sys.argv[3], 'w').write('\n'.join(kept) + '\n')
print('kept', len(kept))
