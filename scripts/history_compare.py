# -*- coding: utf-8 -*-
"""
Author: sacanlee

Model vs. the five historical Kondratiev waves: comparison and evaluation.

Data sources:
- Long-wave periodisation: book (2022) Table 5.1 (idealised long waves) and
  Appendix 8.A (phases of waves 4 and 5: 1948-1966/1967-1982, 1983-2006/2007-2021).
- Population: Our World in Data global population series (1780-2023, based on
  Maddison/HYDE + UN WPP), https://ourworldindata.org/grapher/population.csv;
  recent segment cross-validated with the World Bank API (WLD, SP.POP.GROW).
- n = CAGR of the working-age population (15-64) over each wave's span.

Usage: python history_compare.py
Output: ../data/history_results.json + terminal table
"""
import numpy as np, json, importlib.util, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data')
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location('lwm', os.path.join(HERE, 'long_wave_model.py'))
lwm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lwm)

# ---------- Historical long-wave periodisation (book Table 5.1 + App. 8.A) ----------
# The model input n uses the working-age population (15-64) growth rate:
#   after 1950: UN WPP age-structured data (OWID population-by-age-group.csv,
#               15-64 = Ages 25-64 + Ages 15-24);
#   1940-1950: share linearly interpolated (1940: 59.0%, anchored to the true 1950: 60.1%);
#   1790-1940: shares imputed (1800: 56.0% -> 1900: 57.5% -> 1940: 59.0%, based on
#              19th-century demography with ~40% of the world population aged 0-14;
#              imputation error ±0.03pp).
#   World Bank cross-check: SP.POP.1564.TO (1960+) agrees with OWID (<0.5%);
#   SL.TLF.TOTL.IN (1990+) shows falling labour-force participation (see docs).
# n_pop = original (total population CAGR), n = working-age CAGR; n_rise/n_fall
# likewise (per segment, reference only).
WAVES = [
    dict(name="Wave 1 Industrial", t0=1790, t_peak=1815, t1=1845, T=55, rise=25, fall=30,
         n=0.00537, n_pop=0.00512, n_rise=0.00470, n_fall=0.00594),
    dict(name="Wave 2 Victorian", t0=1845, t_peak=1873, t1=1896, T=51, rise=28, fall=23,
         n=0.00493, n_pop=0.00467, n_rise=0.00334, n_fall=0.00686),
    dict(name="Wave 3 Belle Epoque", t0=1896, t_peak=1920, t1=1940, T=44, rise=24, fall=20,
         n=0.00907, n_pop=0.00845, n_rise=0.00812, n_fall=0.01024),
    dict(name="Wave 4 Golden Age", t0=1940, t_peak=1966, t1=1982, T=42, rise=26, fall=16,
         n=0.01691, n_pop=0.01679, n_rise=0.01622, n_fall=0.01936),
    dict(name="Wave 5 Information", t0=1982, t_peak=2007, t1=2022, T=40, rise=26, fall=15,
         n=0.01623, n_pop=0.01393, n_rise=0.01897, n_fall=0.01169),
]

def model_for_wave(n):
    m = lwm.Model(n=n)
    eq = m.equilibrium_closed_form()
    J, w = m.eig_analysis(eq)
    wc = sorted(w, key=lambda z: -abs(z.imag))
    dom = wc[0]
    T_pred = 2*np.pi/dom.imag
    # 5% perturbation -> linear rise/fall + nonlinear (6T window) rise/fall
    t, sol = m.linear_sim(0.05*eq, J, T_years=300)
    seg = lwm.cycle_segments(t, sol, pick_index=0)
    T0 = T_pred
    t_nl, sol_nl = m.simulate(eq + 0.05*eq, T_years=min(300, 6*T0))
    seg_r = lwm.cycle_segments(t_nl, sol_nl, pick_index=0)
    seg_sc = lwm.cycle_segments(t_nl, sol_nl, pick_index=2)
    return dict(n=n, eq=[round(float(v), 4) for v in eq],
                lam=f"{dom.real:+.4f}{dom.imag:+.4f}i", T_pred=round(T_pred, 1),
                lin=dict(T=round(seg[0],1), rise=round(seg[1],1), fall=round(seg[2],1)) if seg[0]==seg[0] else None,
                nl_r=dict(T=round(seg_r[0],1), rise=round(seg_r[1],1), fall=round(seg_r[2],1)) if seg_r[0]==seg_r[0] else None,
                nl_sc=dict(T=round(seg_sc[0],1), rise=round(seg_sc[1],1), fall=round(seg_sc[2],1)) if seg_sc[0]==seg_sc[0] else None)

def main():
    out = []
    print(f"{'Wave':<18s}{'n':>7s}{'Actual T':>9s}{'Pred T (lin)':>13s}{'Bias%':>7s}"
          f"{'Lin rise/fall':>15s}{'NL rise/fall (sC)':>19s}{'Actual rise/fall':>17s}")
    for w in WAVES:
        r = model_for_wave(w['n'])
        bias = (r['T_pred'] - w['T'])/w['T']*100
        lin = r['lin']; nl = r['nl_sc']
        lins = f"{lin['rise']}/{lin['fall']}" if lin else "none"
        nls = f"{nl['rise']}({100*nl['rise']/(nl['rise']+nl['fall']):.0f}%)/{nl['fall']}" if nl else "none"
        print(f"{w['name']:<18s}{w['n']*100:6.2f}%{w['T']:6d}{r['T_pred']:10.1f}y{bias:7.0f}%"
              f"{lins:>15s}{nls:>19s}{f'{w['rise']}/{w['fall']}':>17s}")
        out.append(dict(wave=w['name'], period=f"{w['t0']}-{w['t1']}",
                        n_input=w['n'],
                        actual_T=w['T'], actual_rise=w['rise'], actual_fall=w['fall'],
                        **r))
    os.makedirs(DATA_DIR, exist_ok=True)
    with open(os.path.join(DATA_DIR, 'history_results.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\nSaved ../data/history_results.json")

if __name__ == '__main__':
    main()
