# -*- coding: utf-8 -*-
"""
Author: sacanlee

Sub-Saharan Africa (SSA) fertility scenario: TFR falls linearly to 2.0 over
2025-2045 (below 2.1 within 20 years), projected change in the world's
working-age population (15-64) growth after 2040, and re-estimation of the
sixth/seventh long waves.

Data: UN WPP 2024 (medium), via OWID:
  fertility-rate-with-projections.csv        (country TFR estimates + projections)
  population-by-age-group-with-projections.csv (country 15-64 age composition)
SSA per UN M49 definition (Eastern 17 + Middle 9 + Southern 5 + Western 17
= 48 economies) aggregated by country.

Method (simplified but transparent):
  - SSA annual births b(c) ~ Under-5s(c+2)/0.95 (0-4 survival approximation, constant first order)
  - TFR(SSA) = population-weighted average of country TFRs;
    scenario path TFRs(t) = linear from 2025 to 2.0 in 2045, then back to the UN path (q=0)
  - q(c) = 1 - TFRs(c)/TFR_UN(c): reduction ratio per birth cohort (2025-2045)
  - World 15-64(Y) correction: reduction = sum_{c: c+15<=Y<=c+64} q(c)*b(c)*0.93
    (15-64 survival approximation)
  - n_wa from the corrected series' annual growth; wave-stage average n drives
    the long-wave model

Run: python ssa_fertility_scenario.py --data-dir C:/tmp
Output: ../data/ssa_scenario_results.json + terminal tables
"""
import csv, json, os, importlib.util, sys, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data')

SSA_CODES = "BDI COM DJI ERI ETH KEN MDG MWI MUS MOZ RWA SOM SSD TZA UGA ZMB ZWE " \
            "AGO CMR CAF TCD COD COG GNQ GAB STP " \
            "BWA LSO NAM SWZ ZAF " \
            "BEN BFA CPV CIV GMB GHA GIN GNB LBR MLI MRT NER NGA SEN SLE TGO"
SSA_CODES = set(SSA_CODES.split())

def load_age(base):
    rows = csv.DictReader(open(os.path.join(base, 'owan_pop_age_proj.csv'), encoding='utf-8'))
    d = []
    for r in rows:
        if not r.get('Ages 25-64 (Projected)'):
            continue
        d.append(dict(code=r['Code'], y=int(r['Year']),
                      tot=int(float(r['Total (Projected)'] or 0)),
                      wa=int(float(r['Ages 25-64 (Projected)'])) +
                            int(float(r['Under-25s (Projected)'])) -
                            int(float(r['Under-15s (Projected)'])),
                      u5=int(float(r['Under-5s (Projected)'] or 0)),
                      u15=int(float(r['Under-15s (Projected)'] or 0)),
                      u25=int(float(r['Under-25s (Projected)'] or 0))))
    return d

def load_tfr(base):
    rows = csv.DictReader(open(os.path.join(base, 'owo_tfr_proj.csv'), encoding='utf-8'))
    out = {}
    for r in rows:
        y = int(r['Year'])
        v = r['Fertility rate (projections) (Projected)'] or r['Fertility rate (estimates)']
        if v:
            try:
                out[(r['Code'], y)] = float(v)
            except ValueError:
                pass
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data-dir', default=r'C:\tmp',
                    help='directory containing the OWID CSVs (owan_pop_age_proj.csv, owo_tfr_proj.csv)')
    args = ap.parse_args()
    base = args.data_dir

    age = load_age(base)
    tfr = load_tfr(base)
    # ---- SSA aggregation ----
    yrs = list(range(2020, 2101))
    tot = {y: 0 for y in yrs}; u5 = {y: 0 for y in yrs}
    wa_sum = {y: 0 for y in yrs}
    tfr_w = {y: (0.0, 0.0) for y in yrs}
    for r in age:
        if r['code'] in SSA_CODES and r['y'] in tot:
            tot[r['y']] += r['tot']; u5[r['y']] += r['u5']; wa_sum[r['y']] += r['wa']
            if (r['code'], r['y']) in tfr:
                w, wt = tfr_w[r['y']]
                tfr_w[r['y']] = (w + tfr[(r['code'], r['y'])] * r['tot'], wt + r['tot'])
    tfr_un = {y: (tfr_w[y][0]/tfr_w[y][1] if tfr_w[y][1] else np.nan) for y in yrs}
    print('SSA aggregation check: 2025 pop %.2f bn, weighted TFR %.2f | 2050: %.2f bn, TFR %.2f' % (
        tot[2025]/1e8, tfr_un[2025], tot[2050]/1e8, tfr_un[2050]))
    # ---- Scenario TFR ----
    t0, t1, target = 2025, 2045, 2.0
    tfr_sc = {}
    for y in yrs:
        if y <= t0: tfr_sc[y] = tfr_un[y]
        elif y <= t1: tfr_sc[y] = tfr_un[t0] + (target - tfr_un[t0]) * (y - t0)/(t1 - t0)
        else: tfr_sc[y] = tfr_un[y]
    q = {y: max(0.0, 1 - tfr_sc[y]/tfr_un[y]) for y in yrs if tfr_un[y] not in (None, 0, np.nan)}
    def births(c):  # SSA births in year c (backcast from 0-4 population / 0.95)
        return u5.get(c + 2, 0) / 0.95
    # ---- World 15-64 correction ----
    world = {}
    for r in age:
        if r['code'] == 'OWID_WRL':
            world[r['y']] = r['wa']
    reduced_birth_years = [c for c in range(t0, t1 + 1)]
    delta = {}
    SURV = 0.93
    for y in range(2035, 2101):
        d = 0.0
        for c in reduced_birth_years:
            if c + 15 <= y <= c + 64:
                d += q.get(c, 0) * births(c) * SURV
        delta[y] = d
    wa_cf = {y: world[y] - delta.get(y, 0) for y in world}
    def grow(series, a, b):
        return (series[b]/series[a])**(1/(b-a)) - 1
    print('\nWorld 15-64 (UN medium vs SSA scenario): 2040: %.2fbn vs %.2fbn | 2070: %.2fbn vs %.2fbn | 2100: %.2fbn vs %.2fbn' % (
        world[2040]/1e8, wa_cf[2040]/1e8, world[2070]/1e8, wa_cf[2070]/1e8,
        world[2100]/1e8, wa_cf[2100]/1e8))
    print('Cumulative birth gap (2025-2045 scenario): %.1f bn, effective after ~2050 => max 15-64 reduction ~%.1f bn' % (
        sum(q.get(c, 0)*births(c) for c in reduced_birth_years)/1e8,
        max(delta.values())/1e8))
    print('\nWindow n_wa: UN-medium vs SSA-scenario')
    for (a, b) in [(2030, 2040), (2030, 2045), (2045, 2070), (2070, 2085), (2085, 2100), (2040, 2070)]:
        n0 = grow(world, a, b); n1 = grow(wa_cf, a, b)
        print(f'  {a}-{b}: {100*n0:+.3f}% -> {100*n1:+.3f}%  (delta {100*(n1-n0):+.3f}pp)')
    # ---- Model re-estimation ----
    sys.path.insert(0, HERE)
    spec = importlib.util.spec_from_file_location('lwm', os.path.join(HERE, 'long_wave_model.py'))
    lwm = importlib.util.module_from_spec(spec); spec.loader.exec_module(lwm)
    def model(n):
        m = lwm.Model(n=n)
        eq = m.equilibrium_closed_form()
        J, w = m.eig_analysis(eq)
        dom = sorted(w, key=lambda z: -abs(z.imag))[0]
        T = 2*np.pi/dom.imag
        t_nl, sol_nl = m.simulate(eq + 0.05*eq, T_years=min(300, 6*T))
        s = lwm.cycle_segments(t_nl, sol_nl, pick_index=2)
        p = s[1]/(s[1]+s[2])
        return T, p * s[0], T * p
    print('\n=== Wave 6 (from 2030) re-estimation ===')
    n_rise_cf = grow(wa_cf, 2030, 2045); n_fall_cf = grow(wa_cf, 2045, 2070)
    T1, rise_nl, _ = model(n_rise_cf); T2, _, _ = model(n_fall_cf)
    print(f'  n(rise 2030-45)={100*n_rise_cf:+.3f}%  n(fall 2045-70)={100*n_fall_cf:+.3f}%')
    print(f'  T(rise)={T1:.1f} nonlinear sC rise-segment (measure)~={rise_nl:.1f}y')
    # Main measure: eigenvalue period + nonlinear rise share
    m1 = lwm.Model(n=n_rise_cf); eq1 = m1.equilibrium_closed_form()
    J1, w1 = m1.eig_analysis(eq1); dom1 = sorted(w1, key=lambda z: -abs(z.imag))[0]
    T1 = 2*np.pi/dom1.imag
    t1, s1 = m1.simulate(eq1+0.05*eq1, T_years=min(300, 6*T1))
    seg1 = lwm.cycle_segments(t1, s1, pick_index=2)
    p1 = seg1[1]/(seg1[1]+seg1[2])
    m2 = lwm.Model(n=n_fall_cf); eq2 = m2.equilibrium_closed_form()
    J2, w2 = m2.eig_analysis(eq2); dom2 = sorted(w2, key=lambda z: -abs(z.imag))[0]
    T2 = 2*np.pi/dom2.imag
    t2, s2 = m2.simulate(eq2+0.05*eq2, T_years=min(300, 6*T2))
    seg2 = lwm.cycle_segments(t2, s2, pick_index=2)
    p2 = seg2[1]/(seg2[1]+seg2[2])
    rise = p1*T1; fall = (1-p2)*T2
    print(f'  Wave 6: rise {rise:.1f}y | peak ~{2030+rise:.0f} | fall {fall:.1f}y | total {rise+fall:.1f}y | trough ~{2030+rise+fall:.0f}')
    print('\n=== Wave 7 (from 2070) re-estimation ===')
    n_r7 = grow(wa_cf, 2070, 2085); n_f7 = grow(wa_cf, 2085, 2100)
    m1 = lwm.Model(n=n_r7); eq1 = m1.equilibrium_closed_form()
    J1, w1 = m1.eig_analysis(eq1); dom1 = sorted(w1, key=lambda z: -abs(z.imag))[0]
    T1 = 2*np.pi/dom1.imag
    t1, s1 = m1.simulate(eq1+0.05*eq1, T_years=min(300, 6*T1))
    seg1 = lwm.cycle_segments(t1, s1, pick_index=2)
    p1 = seg1[1]/(seg1[1]+seg1[2])
    m2 = lwm.Model(n=n_f7); eq2 = m2.equilibrium_closed_form()
    J2, w2 = m2.eig_analysis(eq2); dom2 = sorted(w2, key=lambda z: -abs(z.imag))[0]
    T2 = 2*np.pi/dom2.imag
    t2, s2 = m2.simulate(eq2+0.05*eq2, T_years=min(300, 6*T2))
    seg2 = lwm.cycle_segments(t2, s2, pick_index=2)
    p2 = seg2[1]/(seg2[1]+seg2[2])
    rise = p1*T1; fall = (1-p2)*T2
    print(f'  n(rise 2070-85)={100*n_r7:+.3f}%  n(fall 2085-2100)={100*n_f7:+.3f}%')
    print(f'  Wave 7: rise {rise:.1f}y | peak ~{2070+rise:.0f} | fall {fall:.1f}y | total {rise+fall:.1f}y | trough ~{2070+rise+fall:.0f}')
    os.makedirs(DATA_DIR, exist_ok=True)
    json.dump(dict(ssa_check=dict(pop2025=tot[2025]/1e8, tfr2025=tfr_un[2025],
                                  pop2050=tot[2050]/1e8, tfr2050=tfr_un[2050]),
                   reduced_births_gap=sum(q.get(c, 0)*births(c) for c in reduced_birth_years)/1e8,
                   max_delta_wa=max(delta.values())/1e8,
                   wa_cf={k: v/1e8 for k, v in wa_cf.items()}),
              open(os.path.join(DATA_DIR, 'ssa_scenario_results.json'), 'w'), ensure_ascii=False, indent=2)

if __name__ == '__main__':
    main()
