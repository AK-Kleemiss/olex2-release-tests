# stage_cod.py <split-file|all> <n> <out-ids> [seed] [skip-file...]: COD mirror entries -> <OUT>/cod_<id>
# COD=<mirror root> (cif/ and hkl/ trees), OUT=<samples root>; reflections from the separate hkl, else from the
# CIF's _shelx_hkl_file block or its own _refln loop (21 Sep 2026)
import random, re, sys, os, shutil, subprocess
random.seed(int(sys.argv[4]) if len(sys.argv) > 4 else 17)
COD = os.environ.get('COD', 'C:/Users/florian/cctbx-work/COD')
OUT = os.environ.get('OUT', 'D:/devel/olex2-samples')
TOOLS = os.path.dirname(os.path.abspath(__file__))
def path(kind, cid): return '%s/%s/%s/%s/%s/%s.%s' % (COD, kind, cid[0], cid[1:3], cid[3:5], cid, kind)
if sys.argv[1] == 'all':
  ids = sorted(set(f[:-4] for k in ('hkl', 'cif') for r, d, fs in os.walk(COD + '/' + k) for f in fs if f[-4:] == '.' + k))
else:
  ids = [l.strip() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
skip = set(l.split()[0].replace('cod_', '') for f in sys.argv[5:] for l in open(f) if l.strip())
random.shuffle(ids)
want, kept = int(sys.argv[2]), []
def embedded(cif):
  m = re.search(r'(?m)^_shelx_hkl_file\s*\n;\s*\n(.*?)^;', cif, re.S)
  if m: return m.group(1), False
  if '_refln_index_h' in cif and ('_refln_F_squared_meas' in cif or '_refln_F_meas' in cif): return cif, True
  return None, False
for cid in ids:
  if len(kept) >= want: break
  if cid in skip or not os.path.isfile(path('cif', cid)): continue
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
  if os.path.isfile(path('hkl', cid)):
    hkl, loop = open(path('hkl', cid), encoding='utf-8', errors='replace').read(), True
    if '_refln_F_squared_meas' not in hkl or '_refln_F_squared_sigma' not in hkl: continue
  else:
    hkl, loop = embedded(cif)
    if hkl is None: continue
  d = '%s/cod_%s' % (OUT, cid)
  os.makedirs(d, exist_ok=True)
  # a refln loop is the fcf list, already in the RES basis: an HKLF matrix would transform it twice;
  # a _shelx_hkl_file block is the raw file the embedded HKLF line belongs to
  if loop: cif = re.sub(r'(?m)^HKLF\s+(\d+).*$', r'HKLF ', cif)
  open('%s/%s.cif' % (d, cid), 'w', encoding='utf-8').write(cif)
  open('%s/%s.hkl' % (d, cid), 'w', encoding='utf-8').write(hkl)
  if loop:
    r = subprocess.run([sys.executable, os.path.join(TOOLS, 'cod2hkl.py'), '%s/%s.hkl' % (d, cid)], capture_output=True, text=True)
    if r.returncode: shutil.rmtree(d); print('drop', cid, r.stderr.strip()[-80:]); continue
  kept.append(cid); print(cid, len(nonh), m.group(1), 'loop' if loop else 'shelx', flush=True)
open(sys.argv[3], 'w').write('\n'.join(kept) + '\n')
print('kept', len(kept))
