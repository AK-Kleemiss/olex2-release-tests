import sys
root = sys.argv[1]


def patch(path, start, end, new):
  s = open(path, newline='').read()
  nl = '\r\n' if '\r\n' in s else '\n'
  a, b = s.index(start), s.index(end)
  assert 0 < a < b
  s = s[:a] + new.replace('\n', nl) + s[b:]
  open(path, 'w', newline='').write(s)


ADAPTER = '''  def refinedElectronCounts(self, xs):
    """ An electron count per atom of the refined model.

    Fo and Fc maps on one grid, both phased by the model and integrated in the
    same sphere about each site: the ratio times the modelled Z is what the
    data say sits there, the thermal smearing and series termination
    cancelling per site. Divided by its median over the better-behaved half
    of the atoms so the scale, the missing hydrogens and the noise peaks do
    not move every atom together. None where the model has no
    element (a Q peak).
    """
    from cctbx import maptbx, miller
    from cctbx.array_family import flex
    from cctbx.eltbx import tiny_pse
    from smtbx.ab_initio import element_assignment
    f_obs = self.reflections.f_sq_obs_merged.average_bijvoet_mates().f_sq_as_f()
    f_calc = f_obs.structure_factors_from_scatterers(xray_structure=xs).f_calc()
    fc = flex.abs(f_calc.data())
    k = flex.sum(f_obs.data()*fc)/flex.sum(fc*fc)
    gridding = maptbx.crystal_gridding(
      unit_cell=xs.unit_cell(), space_group_info=xs.space_group_info(),
      d_min=f_obs.d_min(), resolution_factor=1/3.,
      symmetry_flags=maptbx.use_space_group_symmetry)
    def integrate(data):
      m = miller.fft_map(gridding, miller.array(miller_set=f_calc, data=data))
      m.apply_volume_scaling()
      return element_assignment.integrated_densities(m, xs.sites_frac())
    eo = integrate(flex.polar(f_obs.data()/k, flex.arg(f_calc.data())))
    ec = integrate(f_calc.data())
    ratio = [o/c if c > 0 else None for o, c in zip(eo, ec)]
    # ponytail: the median over the low-U half only; noise peaks (high U, no
    # density) sit low and would drag the median under the real atoms
    u = list(xs.extract_u_iso_or_u_equiv())
    u_mid = sorted(u)[len(u)//2] if u else 0.0
    good = sorted(r for r, ui in zip(ratio, u) if r is not None and ui <= u_mid)
    norm = good[len(good)//2] if good else 1.0
    # ponytail: a model typed a Z too light everywhere is self-consistent
    # under the median (11 O read as C, 2211782); the formula counts, where
    # the user gave real ones, pin the middle half of the Z distribution
    have, want = [], self.expectedZs()
    for s in xs.scatterers():
      try:
        have.append(tiny_pse.table(s.scattering_type.strip().capitalize()).atomic_number())
      except (RuntimeError, ValueError):
        pass
    def mid(v):
      v = sorted(v)
      v = v[len(v)//4:len(v) - len(v)//4] or v
      return sum(v)/float(len(v)) if v else 0.0
    # ponytail: rank-based, so it only holds when the model and the formula
    # count about the same atoms (2223999: 9 modelled against 21 declared
    # read Mo and K as O); the window is the knob
    if want and 0.75 <= len(want)/float(len(have) or 1) <= 1.33 and mid(have) > 0:
      print("Formula anchor: model mid-Z %.1f, formula %.1f, scale x%.2f"
            % (mid(have), mid(want), mid(have)/mid(want)))
      norm *= mid(have)/mid(want)
    out = []
    for s, r in zip(xs.scatterers(), ratio):
      try:
        z = tiny_pse.table(s.scattering_type.strip().capitalize()).atomic_number()
      except (RuntimeError, ValueError):
        z = 0
      out.append(z*r/norm if r is not None and z > 0 else None)
    return out

  def expectedZs(self):
    """ One Z per non-H atom of the declared formula; empty when the counts
    look qualitative (one of each), which SHELXT users type. """
    from cctbx.eltbx import tiny_pse
    out, ns = [], []
    try:
      for part in str(olx.xf.GetFormula('list')).split(','):
        e, n = part.split(':')
        e = e.strip().capitalize()
        if e not in ("H", "D"):
          ns.append(int(round(float(n))))
          out += [tiny_pse.table(e).atomic_number()]*ns[-1]
    except Exception:
      return []
    return out if len(ns) > 1 and max(ns) > 1 or len(ns) == 1 else []

  def completeModel(self):
    """ Add what the difference map of the refined model says is missing.

    smtbx.ab_initio.model_completion: rounds of Fo-Fc on the refined model,
    peaks over its sigma cut at bond distance to the model become carbon
    atoms, then a refine and a prune of what it added. Returns the number
    of atoms posted; never raises.
    """
    import re
    try:
      from smtbx.ab_initio import model_completion
      xs = self.xray_structure()
      f_obs = self.reflections.f_sq_obs_merged.average_bijvoet_mates().f_sq_as_f()
      scat = [s for s in xs.scatterers()
              if s.scattering_type.strip().capitalize() not in ("Q", "H", "D")]
      sites = [s.site for s in scat]
      calls = [s.scattering_type.strip().capitalize() for s in scat]
      new_sites, new_calls, n = model_completion.complete(f_obs, xs, sites, calls)
      if n <= 0:
        print("Difference map completion: nothing to add")
        return 0
      used = [int(m.group(1)) for m in
              (re.match(r"C(\d+)$", str(s.label)) for s in xs.scatterers()) if m]
      next_no = max(used or [0]) + 1
      posted = []
      for k, site in enumerate(list(new_sites)[len(sites):]):
        name = self.post_single_peak(site, 0, element="C", number=next_no + k)
        if name:
          posted.append(name)
      print("Difference map completion added %d atom(s): %s"
            % (len(posted), " ".join(posted)))
      return len(posted)
    except Exception as e:
      import traceback
      print("Difference map completion did not run (%s: %s)"
            % (type(e).__name__, e))
      if OV.IsDebugging():
        traceback.print_exc()
      return 0

  def printDoubt(self):
    """ The atoms whose type the last re-typing did not settle: every other
    allowed element that kept over a tenth of the arbitration weight. """
    rows = getattr(type(self), "doubt", None)
    if not rows:
      return
    print("Types with a runner-up over 10 %:")
    for name, alt in rows:
      print("  %s: %s" % (name, " | ".join("%s (%.0f %%)" % (e, 100*p)
                                            for e, p in alt)))

  def reassignAfterCleanup(self):
    """ Re-type the refined model from what the data say sits at each site.

    Density evidence is taken fresh from the refined model rather than from
    the solution's peaks: an atom typed too light shows the surplus in the Fo
    map, one typed too heavy a deficit, and a correctly typed atom reads its
    own Z. An atom that reads its own Z is left alone, whatever the geometry
    classifier thinks; the others are arbitrated between a Gaussian in Z
    about the estimate and the classifier's posterior (floored so that a
    class the classifier never saw still gets the density's vote), over
    the declared composition. Returns {label: (old, new)} for the labels changed; never
    raises.
    """
    import math
    from cctbx.array_family import flex
    from cctbx.eltbx import tiny_pse

    allowed = set(_SOLUTION_EVIDENCE.get("allowed") or self.expectedElements())
    # a hydrogen is never what a peak with too little density is
    allowed -= set(("H", "D"))
    try:
      xs = self.xray_structure()
      z_of = [(e, tiny_pse.table(e).atomic_number()) for e in sorted(allowed)]
      # ponytail: a mistyped atom mis-reads its correction (1.6x over a Z away,
      # under once its U has collapsed), so it is re-typed to the nearest Z and
      # read again until it settles; the last two reads averaged, which puts a
      # site hopping between two neighbours in between (four sites measured)
      work, reads = xs.deep_copy_scatterers(), []
      for it in range(4):
        reads.append(self.refinedElectronCounts(work))
        moved = False
        for s, z in zip(work.scatterers(), reads[-1]):
          e = min(z_of, key=lambda ez: abs(ez[1] - z))[0] if z else None
          if e and e != s.scattering_type.strip().capitalize():
            s.scattering_type, moved = e, True
        if not moved:
          break
        work.discard_scattering_type_registry()
      names, sites, z_est, us = [], flex.vec3_double(), [], []
      u_all = xs.extract_u_iso_or_u_equiv()
      u_med = sorted(u_all)[len(u_all)//2] if len(u_all) else 0.0
      for s, z, w, u in zip(xs.scatterers(), reads[-1],
                            reads[-2 if len(reads) > 1 else -1], u_all):
        z = (z + w)/2 if z and w else z
        if z is None:
          continue
        names.append(str(s.label))
        sites.append(s.site)
        z_est.append(z)
        us.append(u)
      if not names:
        return {}
      got = self.geometryProposals(xs.unit_cell(), xs.space_group(), sites)
      proposals = got[0] if got else []
      current = dict((str(s.label), s.scattering_type.strip().capitalize())
                     for s in xs.scatterers())
      changed, doubt = {}, []
      for j, (name, z, u) in enumerate(zip(names, z_est, us)):
        now = current[name]
        z_now = tiny_pse.table(now).atomic_number()
        heavier = [ez for ez in z_of if ez[1] > z_now]
        top = dict(proposals[j]) if j < len(proposals) else {}
        sigma = max(0.35, 0.04*z)
        score = dict((e, math.exp(-0.5*((z - ez)/sigma)**2)
                      * max(top.get(e, 0.0), 0.05)**0.5) for e, ez in z_of)
        tot = sum(score.values()) or 1.0
        # ponytail: a missing heavy atom reads far under its Z (Mo typed O read
        # 18) or its collapsed U absorbs the surplus and it reads its own Z (Ru
        # typed Cl at a fifth of the median U): either steps to the next
        # heavier element, the next round refines it and reads again
        # ponytail: a model whose median U sits under 0.01 is degenerate (the
        # wrong space group doubles every atom) and every U reads collapsed
        if heavier and (z > 1.5*z_now or 0.01 < u_med and u < 0.3*u_med):
          changed[name] = (now, min(heavier, key=lambda ez: ez[1])[0])
        # ponytail: a read within 0.3 of the own Z is left alone; a true N
        # typed C reads 6.4-6.9 and so do some carbons, so in that band the
        # Gaussian is flat over a Z and the geometry posterior decides
        elif abs(z - z_now) >= max(0.3, 0.03*z):
          best = max(score, key=score.get)
          if best != now:
            changed[name] = (now, best)
        final = changed.get(name, (now, now))[1]
        alt = [(e, p/tot) for e, p in sorted(score.items(), key=lambda ep: -ep[1])
               if e != final and p/tot >= 0.1]
        if alt:
          doubt.append((name, alt))
      type(self).doubt = doubt
      if not changed:
        print("Re-typed %d atoms after cleanup; no label changed" % len(names))
        self.printDoubt()
        return changed
      groups = {}
      for name, (old, new) in changed.items():
        groups.setdefault(new, []).append(name)
      for symbol, group in sorted(groups.items()):
        olex.m("sel %s" % " ".join(group))
        olex.m("name sel %s" % symbol)
        olex.m("sel -u")
      print("Re-typed %d atoms after cleanup; %d label(s) changed: %s"
            % (len(names), len(changed),
               ", ".join("%s %s->%s" % (n, o, w)
                         for n, (o, w) in sorted(changed.items()))))
      return changed
    except Exception as e:
      import traceback
      print("Re-typing after cleanup did not run (%s: %s); the labels are "
            "unchanged" % (type(e).__name__, e))
      if OV.IsDebugging():
        traceback.print_exc()
      return {}

'''
patch(root + '/CctbxLib/cctbx_olex_adapter.py',
      '  def refinedElectronCounts(self, xs):', '  def multiTrialSolution(', ADAPTER)

TIDY = '''  def tidy_solution(self):
    """ Make the solution look like chemistry, then drop what will not refine.

    `compaq -a` gathers the fragments, then rounds of: refine `tidy_cycles`,
    prune every atom whose U exceeds `tidy_uiso_factor` times the median U
    of the model (a noise peak has no density to hold it and its U runs
    away; the median follows the temperature and the data where a fixed
    expectation per element did not), gather the difference peaks, make the
    atoms over `tidy_anis_z` anisotropic. Once a
    round (not the first) prunes nothing the survivors are re-typed from the
    refined model and the fragments gathered again (noise peaks typed carbon steal scale from every real atom
    and read them all a Z too heavy), at most `tidy_passes` times. Stops once
    a round changes nothing. Never lets the tidy-up cost the user the
    solution. With `complete_after_tidy` the difference map is asked once for
    what is missing (`completeModel`) once the types are settled.
    """
    factor = OV.GetParam('snum.solution.tidy_uiso_factor') or 3.0
    cycles = OV.GetParam('snum.solution.tidy_cycles')
    cycles = 4 if cycles is None else int(cycles)
    passes = OV.GetParam('snum.solution.tidy_passes')
    passes = 4 if passes is None else int(passes)
    retype = OV.GetParam('snum.solution.retype_after_tidy')
    retype = retype is None or retype
    # ponytail: off by default; on 45 COD entries the additions were all
    # pruned again, the misses were wrong-SG or noise-displaced, not absent
    complete = bool(OV.GetParam('snum.solution.complete_after_tidy'))
    anis_z = OV.GetParam('snum.solution.tidy_anis_z')
    anis_z = 10 if anis_z is None else int(anis_z)
    from cctbx.eltbx import tiny_pse
    def heavy_z(t):
      try:
        return tiny_pse.table(t.capitalize()).atomic_number()
      except (RuntimeError, ValueError):
        return 0
    try:
      olex.m("compaq -a")
    except Exception as err:
      print("Could not assemble the fragments: %s" % err)
    if cycles <= 0:
      return
    try:
      from cctbx_olex_adapter import OlexCctbxSolve
      # ponytail: prune-only rounds do not count, the budget is for re-typing;
      # the first round never re-types, the noise is always still there
      for p in range(3*passes):
        olex.m("refine %d" % cycles)
        us = []
        for i in range(int(olx.xf.au.GetAtomCount())):
          if olx.xf.au.IsAtomDeleted(i) == 'true' or \\
             str(olx.xf.au.GetAtomType(i)) in ('Q', 'H'):
            continue
          try:
            us.append((float(olx.xf.au.GetAtomUiso(i)), i,
                       str(olx.xf.au.GetAtomName(i)),
                       str(olx.xf.au.GetAtomType(i))))
          except (TypeError, ValueError):
            pass
        us.sort()
        median = us[len(us)//2][0] if us else 0.0
        doomed = [n for u, i, n, t in us if u > factor*max(median, 0.005)]
        # ponytail: the heavy atoms go anisotropic once the noise is gone;
        # an isotropic Pd leaves a residual the light atoms then read
        heavy = sorted(set("$" + t for u, i, n, t in us if anis_z
                           and heavy_z(t) > anis_z))
        if doomed:
          olex.m("kill %s" % " ".join(doomed))
          print("Pruned %d peak(s) whose U exceeded %.1fx the median U %.3f "
                "after %d cycles (%d left): %s"
                % (len(doomed), factor, median, cycles, len(us) - len(doomed),
                   " ".join(doomed[:12]) + (" ..." if len(doomed) > 12 else "")))
        else:
          print("Every atom refined to within %.1fx the median U %.3f; "
                "nothing pruned" % (factor, median))
        olex.m("compaq -q")
        if heavy:
          olex.m("anis %s" % " ".join(heavy))
        if doomed or p == 0:
          continue
        # re-typed only once nothing was pruned: noise peaks typed carbon steal
        # scale from every real atom and read them all a Z too heavy
        if not retype or not OlexCctbxSolve().reassignAfterCleanup():
          # ponytail: one completion once the types are settled; what it adds
          # goes through the same refine, prune and re-typing as the rest
          if not complete or not OlexCctbxSolve().completeModel():
            break
          complete = False
          continue
        passes -= 1
        # ponytail: the prunes leave fragments scattered again; a second
        # assembly after the re-typing joins what the first round could not
        olex.m("compaq -a")
        if passes <= 0:
          olex.m("refine %d" % cycles)
          OlexCctbxSolve().printDoubt()
          break
    except Exception as err:
      import traceback
      print("Post-solution tidy-up failed: %s" % err)
      if OV.IsDebugging():
        traceback.print_exc()

'''
patch(root + '/PyToolLib/method_imp/cctbx.py',
      '  def tidy_solution(self):', '  def show_space_group_suggestions(', TIDY)
print("patched")

NORM = '''    formula = {}
    for element in str(olx.xf.GetFormula('list')).split(','):
      element_type, n = element.split(':')
      formula.setdefault(element_type, float(n))
    if params.amplitude_type == 'E':
      extra.normalisations_for = lambda f: f.amplitude_normalisations(formula)
    elif params.amplitude_type == 'quasi-E':
      def quasi(f):
        try:
          return charge_flipping.amplitude_quasi_normalisations(f)
        except AssertionError:
          # ponytail: the binned means extrapolate under zero on a few hundred
          # reflections (2222909: 311); the formula's Wilson E instead
          print("quasi-E normalisation failed on %d reflections; using E from "
                "the formula" % f.size())
          return f.amplitude_normalisations(formula)
      extra.normalisations_for = quasi
'''
patch(root + '/CctbxLib/cctbx_olex_adapter.py',
      "    if params.amplitude_type == 'E':", chr(10) + "    # Set on every run so a previous run's table", NORM)

patch(root + '/CctbxLib/cctbx_olex_adapter.py',
      "                                           elements=sorted(allowed) or None,",
      "      calls = assigned.assignments",
      """                                           elements=sorted(allowed) or None,
                                           space_group=f_calc.space_group(),
                                           formula_zs=self.expectedZs() or None)
      print("Density scale from %s" % assigned.scale_from)
""")
