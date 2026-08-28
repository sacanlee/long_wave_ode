# -*- coding: utf-8 -*-
"""
Author: sacanlee

2050-2080 global working-age population (15-64) growth: three scenarios x
long-wave model (nonlinear) forecast.

Data and scenarios (2026-08-23):
- Medium: UN WPP 2024 (medium) age-structured projections, OWID chart
  "population-by-age-group-with-projections" (UN WPP 2024).
  15-64 = Ages 25-64 + (Under-25s - Under-15s).
  World 15-64 in 2050: 6.117bn (63.3% of total), peak 2070: 6.286bn,
  2080: 6.205bn (60.3%). 2050-2080 n(15-64) = +0.048%/yr (n_pop = +0.207%/yr,
  share term -0.159pp).
- Optimistic/Pessimistic: UN WPP 2024 High/Low fertility variants (UN published
  totals: High 10.4/12.5bn, Medium 9.7/10.3bn, Low 9.1/8.7bn for 2050/2080).
  n_wa(variant) ~ n_wa(medium) + [n_pop(variant) - n_pop(medium)]
  (assumes the age-structure drift is similar across variants - conservative;
  the Low variant is actually more pessimistic, direction left open).
- Lancet/IHME (Vollset et al. 2024, Lancet GBD fertility): reference scenario
  peaks at ~9.73bn in 2064, ~8.79bn in 2100 - clearly below UN medium (10.3bn
  in 2080), i.e. between UN medium and low -> its reference scenario is roughly
  "pessimistic-leaning-medium" (cross-referenced in the docs).

Run: python forecast_2050_2080.py
Output: ../data/forecast_results.json + terminal table
"""
import numpy as np, importlib.util, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data')
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location('lwm', os.path.join(HERE, 'long_wave_model.py'))
lwm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lwm)

# ---------- Scenario n (annual, decimal) ----------
SCEN = [
    dict(name="Optimistic (UN High fertility variant)", n=+0.00464,
         note="UN High: world total population 10.4/12.5bn in 2050/2080; n_wa ~ +0.048% + (0.617%-0.200%)"),
    dict(name="Medium (UN WPP 2024 Medium)", n=+0.00048,
         note="UN Medium actual 15-64 data: 2050-2080 CAGR = +0.048%/yr; turns negative after the 2070 peak"),
    dict(name="Pessimistic (UN Low fertility variant)", n=-0.00302,
         note="UN Low: 9.1/8.7bn in 2050/2080; n_wa ~ +0.048% + (-0.150%-0.200%); Lancet/IHME reference scenario lies between medium and low"),
]
# Decade sub-segments of the medium case (UN Medium actual data)
DECADE = [("2050-2060", +0.00162), ("2060-2070", +0.00110), ("2070-2080", -0.00129)]

def run(m, label):
    eq = m.equilibrium_closed_form()
    J, w = m.eig_analysis(eq)
    wc = sorted(w, key=lambda z: -abs(z.imag))
    dom = wc[0]
    T = 2*np.pi/dom.imag
    t, sol = m.linear_sim(0.05*eq, J, T_years=300)
    seg = lwm.cycle_segments(t, sol, pick_index=0)
    t_nl, sol_nl = m.simulate(eq + 0.05*eq, T_years=min(300, 6*T))
    seg_r = lwm.cycle_segments(t_nl, sol_nl, pick_index=0)
    seg_sc = lwm.cycle_segments(t_nl, sol_nl, pick_index=2)
    def f(s):
        return (None if s[0] != s[0] else
                dict(T=round(s[0],1), rise=round(s[1],1), fall=round(s[2],1),
                     rise_pct=round(100*s[1]/(s[1]+s[2]),1)))
    out = dict(label=label, n=round(m.p['n'], 5),
               eq=[round(float(v), 4) for v in eq],
               lam=f"{dom.real:+.4f}{dom.imag:+.4f}i", T_eig=round(T, 1),
               lin=f(seg), nl_r=f(seg_r), nl_sC=f(seg_sc))
    print(f"{label}  n={out['n']*100:+.3f}%/yr")
    print(f"   Equilibrium r*={out['eq'][0]:.4f} sV*={out['eq'][1]:.4f} sC*={out['eq'][2]:.4f} "
          f"delta*={out['eq'][3]:.4f} tau*={out['eq'][4]:.4f}")
    print(f"   Dominant eigenvalue {out['lam']}  -> T (eigenvalue) = {out['T_eig']} years")
    if out['lin']:
        g = out['lin']
        print(f"   Linear: peak-to-peak {g['T']} y | rise {g['rise']}/fall {g['fall']} (rise share {g['rise_pct']}%)")
    for key, nm in (('nl_sC', 'Nonlinear sC'), ('nl_r', 'Nonlinear r')):
        g = out[key]
        if g:
            print(f"   {nm}: peak-to-peak {g['T']} y | rise {g['rise']}({g['rise_pct']}%) / fall {g['fall']}")
    return out

def main():
    res = []
    for s in SCEN:
        m = lwm.Model(n=s['n'])
        r = run(m, s['name'] + f"  [{s['note']}]")
        res.append(r)
    print("\nMedium case by decade (UN Medium actual values):")
    for nm, n in DECADE:
        m = lwm.Model(n=n)
        r = run(m, f"{nm}" + " (decade average)")
        res.append({**r, 'decade': nm})
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(os.path.join(DATA_DIR, 'forecast_results.json'), 'w', encoding='utf-8') as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("\nSaved ../data/forecast_results.json")

if __name__ == '__main__':
    main()
