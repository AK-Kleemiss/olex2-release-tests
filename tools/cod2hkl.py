# cod2hkl.py <file.hkl>...: COD fcf-style reflection file -> SHELX HKLF 4 in place (original kept as .fcf)
import sys, os, shutil
for p in sys.argv[1:]:
  lines = open(p).read().splitlines()
  cols, rows, inloop = [], [], False
  for l in lines:
    t = l.split()
    if not t: continue
    if t[0] == 'loop_':
      if rows: break  # a later loop (embedded CIF) would reset the columns
      inloop, cols = True, []; continue
    if inloop and t[0].startswith('_refln'): cols.append(t[0]); continue
    if inloop and cols and not t[0].startswith('_') and len(t) >= len(cols): rows.append(t)
    elif inloop and cols and t[0].startswith('_'): inloop = False
  h, k, l = (cols.index('_refln_index_' + x) for x in 'hkl')
  sq = '_refln_F_squared_meas' in cols
  f = cols.index('_refln_F_squared_meas' if sq else '_refln_F_meas'); s = cols.index('_refln_F_squared_sigma' if sq else '_refln_F_sigma')
  if not sq:  # an F loop: F^2 and sigma(F^2) = 2 F sigma(F)
    for r in rows: F, S = float(r[f]), float(r[s]); r[f], r[s] = str(F*F), str(2*abs(F)*S)
  shutil.copy(p, p[:-4] + '.fcf')
  # one common scale keeps every value inside the 8-character HKLF 4 field
  scale = 1.0
  while max(abs(float(r[c])) for r in rows for c in (f, s))/scale > 99999.99: scale *= 10
  with open(p, 'w') as o:
    for r in rows:
      o.write('%4d%4d%4d%8.2f%8.2f\n' % (int(r[h]), int(r[k]), int(r[l]), float(r[f])/scale, float(r[s])/scale))
    o.write('   0   0   0    0.00    0.00\n')
  print(p, len(rows), 'reflections')
