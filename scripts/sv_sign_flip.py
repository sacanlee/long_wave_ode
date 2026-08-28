# -*- coding: utf-8 -*-
"""
Author: sacanlee

Test: if the sign of sV in Eq. (8.15) is changed from "-" to "+", are the
internal inconsistencies reported earlier (see Technical_Report_Ch8_Model_Errata.docx,
six findings) "rescued"?

Printed Eq. (8.15):  rdot = -a1*(sC - sV)*r^2 + a2*(delta + tau - sV)*r
Four variants (only the sV signs in (8.15) are changed):
  orig:  -a1(sC - sV)r^2 + a2(delta + tau - sV)r     (as printed)
  A:     -a1(sC + sV)r^2 + a2(delta + tau + sV)r     (both sV signs flipped)
  B:     -a1(sC + sV)r^2 + a2(delta + tau - sV)r     (first bracket only)
  C:     -a1(sC - sV)r^2 + a2(delta + tau + sV)r     (second bracket only)

For each variant:
  F1:  Is the book's claimed equilibrium {0.0947,0.1188,0.8708,0.05,0.045} an equilibrium? (residuals)
  F2:  Does the printed closed form (8.25) reproduce the claimed equilibrium?
  F3:  Is the printed closed form (8.25) correct (vs. the variant's correct closed form)?
  F4:  Is there a positive (economically meaningful) equilibrium?
  F5:  Eigenvalues of the correct Jacobian: stability / period 2*pi/Im (book claims ~50 y, stable node-focus)
  F6:  Nonlinear simulation period (6T window, median)

Output: ../data/sv_flip_results.json + terminal report
"""
import numpy as np
from scipy.integrate import odeint
import sympy as sp
import math, json, os, sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data')

# ---------------- Parameters ----------------
A1, A2, B0, B1, B2, GW, N = 0.04, 0.01, 0.004, 0.005, 0.05, 0.03, 0.015
BOOK_EQ = [0.0947, 0.1188, 0.8708, 0.05, 0.045]
BOOK_MATRIX = np.array([
    [0.0,     -0.0825,  -0.1034,  0.2177,   0.0],
    [0.1263,   0.0,      1.1611,  -1.0,     0.0],
    [-0.00027, 0.00122, -0.0021,  0.00095,  0.00095],
    [0.0,      0.0,      0.0,      0.0,     0.00226],
    [-0.0063,  0.0,     -0.05806, 0.0,     -0.08]])
BOOK_EIGS = [-0.07904, -0.0002+0.0947j, -0.0002-0.0947j, -0.0014+0.0001j, -0.0014-0.0001j]

VARIANTS = {
    'orig': dict(s1=-1.0, s2=-1.0, note='as printed (8.15)'),
    'A':    dict(s1=+1.0, s2=+1.0, note='both sV signs flipped to +'),
    'B':    dict(s1=+1.0, s2=-1.0, note='first bracket only'),
    'C':    dict(s1=-1.0, s2=+1.0, note='second bracket only'),
}

# ---------------- RHS ----------------
def f(y, s1, s2, n=N):
    r, sV, sC, d, tau = y
    rdot = -A1*(sC + s1*sV)*r**2 + A2*(d + tau + s2*sV)*r
    return np.array([rdot,
                     sC*r - d - tau - n,
                     sC*(tau - sV*r),
                     B0 + d*(GW - sC*r),
                     B1*math.erf(d - B2)])

# ---------------- Correct closed-form equilibrium (per-variant analytical) ----------------
def correct_equilibrium(s1, s2, n=N):
    """From 8.16-8.19: delta*=b2, tau*=K-b2-n (K=gw+b0/b2), sC*r*=K, sV*r*=tau*;
    substituting into (8.15): [a1(K+s1*tau) - a2(b2+tau)]*r* = a2*s2*tau."""
    K = GW + B0/B2
    tau = K - B2 - n
    d = B2
    denom = A1*(K + s1*tau) - A2*(B2 + tau)
    r = A2*s2*tau/denom
    sV = tau/r
    sC = K/r
    return np.array([r, sV, sC, d, tau])

# ---------------- Printed closed form (8.25) ----------------
def printed_closed_form(n=N):
    D = B0 - B2*(B2 - GW + n)
    N_ = B2*(B2 + n) - A2*(B0 + B2*(GW - n))
    r = A2*D/(A1*N_)
    sV = (A1/A2)*(B2 + n) - B0/B2 - GW + n
    sC = (B0 + B2*GW)*A1*N_/(A2*B2*D)
    return np.array([r, sV, sC, B2, B0/B2 - B2 + GW - n])

# ---------------- Analytical Jacobian ----------------
_syms = sp.symbols('r sV sC d tau')
_r, _sV, _sC, _d, _t = _syms

def _jac_fun(s1, s2):
    rdot = -A1*(_sC + s1*_sV)*_r**2 + A2*(_d + _t + s2*_sV)*_r
    sVdot = _sC*_r - _d - _t - N
    sCdot = _sC*(_t - _sV*_r)
    ddot = B0 + _d*(GW - _sC*_r)
    tdot = B1*(2/sp.sqrt(sp.pi))*(_d - B2)
    J = sp.Matrix([[sp.diff(fx, v) for v in _syms] for fx in (rdot, sVdot, sCdot, ddot, tdot)])
    return sp.lambdify((_r, _sV, _sC, _d, _t), J, 'numpy')

_JACS = {k: _jac_fun(v['s1'], v['s2']) for k, v in VARIANTS.items()}

def jacobian(variant, x):
    return _JACS[variant](x[0], x[1], x[2], x[3], x[4])

def describe(w):
    out = []
    for lam in sorted(w, key=lambda z: -abs(z.imag)):
        im = abs(lam.imag)
        T = 2*np.pi/im if im > 1e-10 else float('inf')
        out.append(dict(lam=f"{lam.real:+.5f}{lam.imag:+.5f}i",
                        T=(round(T, 1) if T < 1e6 else None)))
    return out

# ---------------- Period measurement ----------------
def cycle_segments(t, y, pick_index=0):
    from scipy.signal import find_peaks
    yv = y[:, pick_index]
    pk, _ = find_peaks(yv, distance=len(t)//80)
    tr, _ = find_peaks(-yv, distance=len(t)//80)
    events = sorted([(int(i), 'pk') for i in pk] + [(int(i), 'tr') for i in tr])
    rises, falls = [], []
    for i in range(1, len(events)):
        dt = t[events[i][0]] - t[events[i-1][0]]
        if events[i][1] == 'pk' and events[i-1][1] == 'tr':
            rises.append(dt)
        elif events[i][1] == 'tr' and events[i-1][1] == 'pk':
            falls.append(dt)
    periods = [ri + fa for ri, fa in zip(rises, falls)]
    def med(x): return float(np.median(x)) if len(x) else float('nan')
    return (med(periods), med(rises), med(falls), len(periods))

# ---------------- Full check of one variant ----------------
def check_variant(name, v):
    s1, s2 = v['s1'], v['s2']
    out = dict(variant=name, note=v['note'], params=dict(a1=A1, a2=A2, b0=B0, b1=B1,
                b2=B2, gw=GW, n=N), book_eq=BOOK_EQ)

    # F1: residuals at the claimed equilibrium
    res = f(BOOK_EQ, s1, s2)
    out['F1_residuals'] = [round(float(x), 6) for x in res]
    out['F1_max_res'] = round(float(max(abs(x) for x in res)), 6)

    # F2: printed closed form (8.25)
    pcf = printed_closed_form()
    out['F2_printed_closed_form'] = [round(float(x), 5) for x in pcf]
    out['F2_matches_book'] = all(abs(a-b) < 0.0005 for a, b in zip(BOOK_EQ, pcf))

    # F3: printed closed form vs the variant's correct closed form
    cq = correct_equilibrium(s1, s2)
    out['F3_correct_equilibrium'] = [round(float(x), 5) for x in cq]
    out['F3_printed_vs_correct'] = [round(float(a-b), 5) for a, b in zip(pcf, cq)]
    out['F3_printed_is_correct'] = all(abs(a-b) < 1e-6 for a, b in zip(pcf, cq))

    # F4: positive equilibrium?
    pos = all(x > 0 for x in cq)
    ssum = cq[1] + cq[2]
    out['F4_positive_eq'] = pos
    out['F4_sC_plus_sV'] = round(float(ssum), 4)
    out['F4_econ_meaningful'] = bool(pos and cq[0] < 1 and ssum <= 1.0001)

    # F5: Jacobian eigenvalues (a) at the variant's correct equilibrium (b) at the claimed equilibrium
    res5 = {}
    for tag, pt in (('correct_eq', cq), ('book_eq', BOOK_EQ)):
        J = jacobian(name, pt)
        w = np.linalg.eigvals(J)
        desc = describe(w)
        dom = next((d0 for d0 in desc if d0['T'] is not None), None)
        stable = all(z.real < 0 for z in w)
        res5[tag] = dict(eigs=desc,
                         dominant=(dom['lam'], dom['T']) if dom else None,
                         stable=stable,
                         period_about_50y=(dom and dom['T'] and 45 < dom['T'] < 55))
    out['F5_jacobian'] = res5

    # F5b: is the printed matrix a numerical approximation of this variant's Jacobian?
    for tag, pt in (('correct_eq', cq), ('book_eq', BOOK_EQ)):
        J = jacobian(name, pt)
        d = np.max(np.abs(J - BOOK_MATRIX))
        res5[tag]['max_abs_diff_vs_printed_matrix'] = round(float(d), 4)
    out['F5_printed_matrix_is_jacobian'] = any(
        res5[t]['max_abs_diff_vs_printed_matrix'] < 0.01 for t in ('correct_eq', 'book_eq'))

    # F6: nonlinear simulation (5% perturbation, 6T window, 300 y cap)
    res6 = {}
    if pos:
        T0 = res5['correct_eq']['dominant'][1] if res5['correct_eq']['dominant'] else 50.0
        win = min(300.0, 6*T0)
        t = np.linspace(0, win, 20001)
        sol = odeint(lambda y, t: f(y, s1, s2), 1.05*cq, t, rtol=1e-9, atol=1e-12)
        for idx, nm in ((0, 'r'), (1, 'sV'), (2, 'sC'), (3, 'delta'), (4, 'tau')):
            P, R, Fa, ncyc = cycle_segments(t, sol, idx)
            if P == P and ncyc >= 1:
                res6[nm] = dict(T=round(P, 1), rise=round(R, 1), fall=round(Fa, 1),
                                rise_pct=round(100*R/(R+Fa), 1), n_cycles=ncyc)
    out['F6_nonlinear'] = res6

    # Printed eigenvalues vs printed matrix (recomputed, baseline self-check)
    wm = np.linalg.eigvals(BOOK_MATRIX)
    out['printed_matrix_eigs_recomputed'] = [f"{z.real:+.5f}{z.imag:+.5f}i" for z in
                                             sorted(wm, key=lambda z: -abs(z.imag))]
    return out

# ---------------- Historical comparison + forecast (flipped model = variant A/C) ----------------
WAVES = [
    dict(name="W1 1790-1845", n=0.00537, T=55), dict(name="W2 1845-1896", n=0.00493, T=51),
    dict(name="W3 1896-1940", n=0.00907, T=44), dict(name="W4 1940-1982", n=0.01691, T=42),
    dict(name="W5 1982-2022", n=0.01623, T=40),
]
SCEN = [
    dict(name="Optimistic (UN High)", n=+0.00464), dict(name="Medium (UN Medium)", n=+0.00048),
    dict(name="Pessimistic (UN Low)", n=-0.00302),
]

def run_model(s1, s2, n):
    K = GW + B0/B2
    tau = K - B2 - n
    d = B2
    denom = A1*(K + s1*tau) - A2*(B2 + tau)
    r = A2*s2*tau/denom
    eq = np.array([r, tau/r, K/r, d, tau])
    vname = ('A' if (s1, s2) == (1.0, 1.0) else 'C' if (s1, s2) == (-1.0, 1.0)
             else 'B' if (s1, s2) == (1.0, -1.0) else 'orig')
    J = jacobian(vname, eq)
    w = np.linalg.eigvals(J)
    dom = sorted(w, key=lambda z: -abs(z.imag))[0]
    T = 2*np.pi/dom.imag
    win = min(300.0, 6*T)
    t = np.linspace(0, win, 20001)
    sol = odeint(lambda y, t: f(y, s1, s2, n), 1.05*eq, t, rtol=1e-9, atol=1e-12)
    seg_r = cycle_segments(t, sol, 0)
    seg_sc = cycle_segments(t, sol, 2)
    def g(s):
        return (None if s[0] != s[0] or s[3] < 1 else
                dict(T=round(s[0], 1), rise=round(s[1], 1), fall=round(s[2], 1),
                     rise_pct=round(100*s[1]/(s[1]+s[2]), 1)))
    return dict(eq=[round(float(x), 4) for x in eq],
                dom_lam=f"{dom.real:+.4f}{dom.imag:+.4f}i", T_eig=round(T, 1),
                stable=bool(all(z.real < 0 for z in w)),
                nl_r=g(seg_r), nl_sC=g(seg_sc))

def main():
    report = {'book_eq': BOOK_EQ, 'printed_matrix': BOOK_MATRIX.tolist(),
              'book_eigs': [str(z) for z in BOOK_EIGS], 'variants': {}}
    print("=" * 100)
    print("Variant check: flipping the sV signs in 8.15 vs the six reported problems")
    print("=" * 100)
    for name, v in VARIANTS.items():
        out = check_variant(name, v)
        report['variants'][name] = {k: val for k, val in out.items()}
        print(f"\n--- Variant {name}: {v['note']} ---")
        print(f"  F1 residuals at claimed equilibrium [r',sV',sC',delta',tau'] = {out['F1_residuals']}  (max={out['F1_max_res']})")
        print(f"  F2 printed closed form (8.25) = {out['F2_printed_closed_form']}  vs book {BOOK_EQ}  match? {out['F2_matches_book']}")
        print(f"  F3 correct closed form (this variant) = {out['F3_correct_equilibrium']}  diff vs printed = {out['F3_printed_vs_correct']}")
        print(f"  F4 positive equilibrium? {out['F4_positive_eq']}  sC*+sV* = {out['F4_sC_plus_sV']}  economically meaningful? {out['F4_econ_meaningful']}")
        r5 = out['F5_jacobian']['correct_eq']
        print(f"  F5 Jacobian at correct equilibrium: dominant {r5['dominant']}  stable? {r5['stable']}  "
              f"period~50y? {r5['period_about_50y']}")
        print(f"     max diff vs printed matrix: correct eq {r5['max_abs_diff_vs_printed_matrix']} | "
              f"claimed eq {out['F5_jacobian']['book_eq']['max_abs_diff_vs_printed_matrix']}")
        r6 = out['F6_nonlinear']
        if r6:
            print(f"  F6 nonlinear period (6T window): " + "  ".join(
                f"{k}={v['T']}y(rise{v['rise_pct']}%)" for k, v in r6.items()))
        else:
            print("  F6 nonlinear: no positive equilibrium, skipped")
    print("\n" + "=" * 100)
    print("Historical five-wave comparison (flipped model, variants A and C; vs printed orig)")
    print("=" * 100)
    hist = {}
    for w in WAVES:
        row = dict(wave=w['name'], actual_T=w['T'])
        for vn in ('orig', 'A', 'C'):
            r = run_model(VARIANTS[vn]['s1'], VARIANTS[vn]['s2'], w['n'])
            row[vn] = dict(T_eig=r['T_eig'], stable=r['stable'], nl_sC=r['nl_sC'],
                           nl_r=r['nl_r'])
        hist[w['name']] = row
        print(f"{w['name']:<16s} actual T={w['T']:3d}y | orig: T_eig={row['orig']['T_eig']:5.1f} "
              f"(sC {row['orig']['nl_sC']}) | A: T_eig={row['A']['T_eig']:5.1f} "
              f"(sC {row['A']['nl_sC']}) | C: T_eig={row['C']['T_eig']:5.1f} (sC {row['C']['nl_sC']})")
    report['history'] = hist
    print("\n" + "=" * 100)
    print("2050-2080 three-scenario forecast (flipped model, variants A and C)")
    print("=" * 100)
    fc = {}
    for s in SCEN:
        row = dict(name=s['name'], n=s['n'])
        for vn in ('orig', 'A', 'C'):
            r = run_model(VARIANTS[vn]['s1'], VARIANTS[vn]['s2'], s['n'])
            row[vn] = dict(eq=r['eq'], T_eig=r['T_eig'], stable=r['stable'],
                           nl_sC=r['nl_sC'], nl_r=r['nl_r'])
        fc[s['name']] = row
        print(f"{s['name']:<24s} n={s['n']*100:+.3f}%/yr | A: T_eig={row['A']['T_eig']:5.1f} "
              f"stable={row['A']['stable']} sC period={row['A']['nl_sC']} | "
              f"C: T_eig={row['C']['T_eig']:5.1f} stable={row['C']['stable']} sC period={row['C']['nl_sC']}")
    report['forecast'] = fc

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(os.path.join(DATA_DIR, 'sv_flip_results.json'), 'w', encoding='utf-8') as fp:
        json.dump(report, fp, ensure_ascii=False, indent=2, default=str)
    print(f"\nSaved ../data/sv_flip_results.json")

if __name__ == '__main__':
    main()
