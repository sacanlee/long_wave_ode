# -*- coding: utf-8 -*-
"""
Author: sacanlee

Class-struggle extension of the long-wave (Kondratiev) model of
Chatzarakis, Tsaliki & Tsoulfidis (2022), "Economic Growth and Long Cycles:
A Classical Political Economy Approach", Routledge, Ch. 8.

Motivation (critique addressed)
-------------------------------
In Section 8.2.1.2 (book p. 212) the authors discuss the "3rd Case: ge < 0",
i.e. a negative growth rate of the rate of surplus value, as the case of
profit-squeeze, and they relate the sign of ge to "relations of production
(income distribution, class struggle, exploitation of labor, etc.)".  When the
model is closed in Sections 8.2.1.3 and 8.3, however, the growth rate of the
rate of surplus value is made a purely mechanical function of the state of the
reserve army, ge =~ -sV (following Marx, Capital I ch. 25 via Glombowski 1983,
ge = e0 - e1*sV, with e0 = 0 and e1 = 1).  The complete five-equation model
(8.15)-(8.19) therefore contains no channel through which an *autonomous*,
"above-normal" wage push of the working class (a class struggle that aims
higher than what the investment in variable capital dictates) could operate.

This script (a) re-derives the book equation with an explicit class-struggle
parameter CS and (b) solves the modified ODE for several CS values, measuring
in particular the impact of CS on (i) the length of the cycle and (ii) the
lengths of the upswing (trough->peak) and downswing (peak->trough).

Formalisation
-------------
The class struggle acts on the growth rate of the rate of surplus value:

    ge = -sV - CS

with CS = 0 reproducing the book.  CS is measured in the same units as ge and
sV (per year).  Positive CS = class struggle above the "normal" level dictated
by the labour market (workers push wages more than the reserve army alone
would warrant); negative CS = struggle below the normal level (a demobilised
working class / a successful capital offensive, e.g. the depression of wages
below the value of labour-power, one of Marx's counteracting causes of the
falling rate of profit, Capital III ch. 14).

Because ge enters the book's model only inside the profit-rate equation
(8.15), the complete CS-extended system is (only Eq. 8.15 changes):

    r'   = -a1*(sC - sV)*r**2 + a2*(delta + tau - sV - CS)*r     (8.15-CS)
    sV'  = sC*r - delta - tau - n                                 (8.16)
    sC'  = sC*(tau - sV*r)                                        (8.17)
    delta' = b0 + delta*(gw - sC*r)                               (8.18)
    tau' = b1*Erf(delta - b2)                                     (8.19)

Economic reading of CS > 0: at every labour-market state the rate of surplus
value grows more slowly (or falls faster) than the reserve-army relation
ge = -sV implies - a profit squeeze that also persists when unemployment is
rising (workers resist wage cuts).  CS < 0: wages lag the labour-market norm
in every phase, so the rate of surplus value grows faster than ge = -sV
implies (repression of wages raises the surplus value at any given state of
the reserve army).

Parameters and reference state (book p. 233): a1=0.04, a2=0.01, b0=0.004,
b1=0.005, b2=0.05, n=0.015, gw=0.03; equilibrium from the book's own closed
form (8.25): {r*, sV*, sC*, delta*, tau*} = {0.1756, 0.1650, 0.6263, 0.05,
0.045}.  Known caveats of the book model (previous reports of this
repository): this closed-form point is not a true fixed point of (8.15)-(8.19)
(no positive fixed point exists for the printed parameters, findings F1-F4)
and the amplitude of the nonlinear solution grows about x3 per cycle around
it (unstable focus), so the model leaves its economically meaningful domain
(sV < 0, the repository convention; the first time a share exceeds the whole
surplus value sC+sV>1 happens even earlier) after roughly 110-140 years.
All comparisons below are therefore made on the *same* initial conditions
(the closed-form point x 1.05, the repository convention) inside the window
t < t_dom, where t_dom = first time sV < 0.

Measurement protocol (per CS value)
-----------------------------------
- Linearised reference: eigenvalues of the analytical Jacobian of the
  CS-extended system at the closed-form point -> nominal period T_lin(CS).
- Nonlinear integration from y0 = 1.05*eq over 6*T_lin years (repository
  convention), then peak/trough detection on each variable:
  * median peak-to-peak period (cycle length);
  * median rise (trough->peak, "upswing") and fall (peak->trough,
    "downswing") duration, and the rise share = rise/(rise+fall);
  completed events are counted only while t < t_dom.
- Profitability regime: mean and end-of-window rate of profit inside t_dom;
  implied net change of ln(e), e = rate of surplus value (post-computed from
  ge = -sV - CS), i.e. the implied movement of the wage share (ln e falls
  when the wage share rises).

Usage:
    python class_struggle_ode.py
    python class_struggle_ode.py --cs -0.1 -0.05 -0.02 -0.01 0 0.01 0.02 0.05 0.1
    python class_struggle_ode.py --no-figs
"""
import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks
import math, json, argparse, os, sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data')
FIG_DIR = os.path.join(REPO_ROOT, 'figures')

# ----------------------------------------------------------------------
# Model parameters (book p. 233)
# ----------------------------------------------------------------------
P = dict(a1=0.04, a2=0.01, b0=0.004, b1=0.005, b2=0.05, gw=0.03, n=0.015)
VARS = ['r', 'sV', 'sC', 'delta', 'tau']


def equilibrium_closed_form(p=P):
    """Equilibrium of the book's closed form (Eq. 8.25) evaluated at the
    given parameters.  CS does not appear in the printed closed form; we keep
    the same reference state for every CS in order to isolate the pure CS
    effect (ceteris paribus the state)."""
    a1, a2, b0, b2, gw, n = [p[k] for k in ('a1', 'a2', 'b0', 'b2', 'gw', 'n')]
    D = b0 - b2*(b2 - gw + n)
    N = b2*(b2 + n) - a2*(b0 + b2*(gw - n))
    r = a2*D/(a1*N)
    sV = (a1/a2)*(b2 + n) - b0/b2 - gw + n
    sC = (b0 + b2*gw)*a1*N/(a2*b2*D)
    d = b2
    tau = b0/b2 - b2 + gw - n
    return np.array([r, sV, sC, d, tau])


def f(y, cs, p=P):
    """Right-hand side of the CS-extended model (8.15-CS)-(8.19)."""
    r, sV, sC, d, tau = y
    a1, a2, b0, b1, b2, gw, n = [p[k] for k in ('a1', 'a2', 'b0', 'b1', 'b2', 'gw', 'n')]
    return np.array([
        -a1*(sC - sV)*r**2 + a2*(d + tau - sV - cs)*r,
        sC*r - d - tau - n,
        sC*(tau - sV*r),
        b0 + d*(gw - sC*r),
        b1*math.erf(d - b2),
    ])


def jacobian_num(y, cs, p=P, h=1e-7):
    """Jacobian of the CS-extended system by central finite differences
    (used only for the reference-state eigenvalue analysis)."""
    nv = len(y)
    out = np.empty((nv, nv))
    for i in range(nv):
        yp = y.copy(); ym = y.copy()
        yp[i] += h; ym[i] -= h
        out[:, i] = (f(yp, cs, p) - f(ym, cs, p)) / (2*h)
    return out


def eig_analysis(cs, eq, p=P):
    """Eigenvalues of the linearisation at the reference state; returns the
    dominant complex pair (largest imaginary part) with its period."""
    J = jacobian_num(eq, cs, p)
    w = np.linalg.eigvals(J)
    w = sorted(w, key=lambda z: -z.imag)
    dom = w[0]
    return dict(lam=str(dom), Re=float(dom.real),
                T_years=(2*np.pi/dom.imag if abs(dom.imag) > 1e-12 else None))


def simulate(cs, y0, T_years, npts=None, p=P):
    npts = npts or max(4001, int(T_years * 80))
    t = np.linspace(0, T_years, npts)
    sol = odeint(lambda y, t: f(y, cs, p), y0, t, rtol=1e-9, atol=1e-12)
    return t, sol


# ----------------------------------------------------------------------
# Cycle measurement
# ----------------------------------------------------------------------
def peak_trough_indices(t, yv, dist_years=10.0):
    dist = max(3, int(dist_years / (t[1] - t[0])))
    pk, _ = find_peaks(yv, distance=dist)
    tr, _ = find_peaks(-yv, distance=dist)
    return pk, tr


def cycle_list(t, sol, idx, t_dom, dist_years=10.0):
    """Completed rise/fall/cycle durations (years) of variable idx from the
    alternating peak/trough events that end before t_dom.  A 'cycle' is the
    distance between two successive peaks (or two successive troughs), formed
    by one rise plus the neighbouring fall."""
    pk, tr = peak_trough_indices(t, sol[:, idx], dist_years)
    ev = sorted([(float(t[i]), 'P') for i in pk] + [(float(t[i]), 'T') for i in tr])
    ris, fal = [], []
    for k in range(1, len(ev)):
        if ev[k][0] > t_dom:
            break                      # events are ordered in time
        if ev[k][1] == 'P' and ev[k-1][1] == 'T':
            ris.append(ev[k][0] - ev[k-1][0])
        elif ev[k][1] == 'T' and ev[k-1][1] == 'P':
            fal.append(ev[k][0] - ev[k-1][0])
    n = min(len(ris), len(fal))        # pairs that share the intermediate event
    return ris[:n], fal[:n], [ris[i] + fal[i] for i in range(n)]


def cycle_stats(t, sol, idx, t_dom, dist_years=10.0):
    """Median period / rise / fall / rise share of variable idx over the
    completed events that end before t_dom.  The individual durations are
    returned as well (cycles/rises/falls lists) for transparency."""
    ris, fal, cycles = cycle_list(t, sol, idx, t_dom, dist_years)
    med = lambda x: float(np.median(x)) if len(x) else float('nan')
    rise_m, fall_m = med(ris), med(fal)
    share = 100*rise_m/(rise_m + fall_m) if (rise_m == rise_m and fall_m == fall_m) else float('nan')
    return dict(n_cycles=len(cycles), T_med=med(cycles), rise_med=rise_m,
                fall_med=fall_m, rise_share_pct=share,
                cycles=[round(c, 2) for c in cycles],
                rises=[round(x, 2) for x in ris],
                falls=[round(x, 2) for x in fal])


def domain_exit_time(t, sol, flag='sV'):
    """First time the solution leaves the economically meaningful domain.
    flag='sV' (repository convention): sV < 0, i.e. absolute disinvestment in
    the wage fund (the reserve army grows without bound).
    flag='shares': also a share above the whole surplus value (sV > 1 or
    sC > 1), the book's premise that sC + sV <= 1."""
    for i in range(1, len(t)):
        sV, sC = sol[i, 1], sol[i, 2]
        if sV < 0 or (flag == 'shares' and (sV > 1.0 or sC > 1.0)):
            return float(t[i])
    return float('inf')


# ----------------------------------------------------------------------
# Scenario runner
# ----------------------------------------------------------------------
def run_scenario(cs, p=P, T_factor=6.0):
    eq = equilibrium_closed_form(p)
    eig = eig_analysis(cs, eq, p)
    y0 = 1.05 * eq
    T_lin = eig['T_years']
    T_run = T_factor * (T_lin if T_lin and T_lin == T_lin else 50.0)
    t, sol = simulate(cs, y0, min(T_run, 400.0))
    # economically meaningful window: the repository convention (sV < 0)
    t_dom = domain_exit_time(t, sol, flag='sV')
    # first time a share exceeds the whole surplus value (book premise sC+sV<=1)
    t_dom_shares = domain_exit_time(t, sol, flag='shares')
    # metrics per variable (events completed inside [0, t_dom))
    var_stats = {}
    for idx, nm in enumerate(VARS):
        var_stats[nm] = cycle_stats(t, sol, idx, t_dom)
    # profitability summary inside the meaningful window
    msk = t <= t_dom
    r_mean = float(sol[msk, 0].mean())
    r_min = float(sol[msk, 0].min())
    r_end = float(sol[msk, 0][-1])
    return dict(cs=cs, t_dom=t_dom, t_dom_shares=t_dom_shares, eig=eig,
                var_stats=var_stats, r_mean=r_mean, r_min=r_min, r_end=r_end,
                t=t, sol=sol, msk=msk)


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cs', type=float, nargs='+',
                    default=[-0.10, -0.05, -0.02, -0.01, -0.005, 0.0, 0.005, 0.01, 0.02, 0.05, 0.10])
    ap.add_argument('--out', default=os.path.join(DATA_DIR, 'class_struggle_results.json'))
    ap.add_argument('--no-figs', action='store_true')
    args = ap.parse_args()

    eq = equilibrium_closed_form()
    print('=' * 100)
    print('Class-struggle (CS) extension of the long-wave model (8.15-CS)-(8.19)')
    print('Reference state = closed-form equilibrium (8.25):',
          ', '.join(f'{nm}*={v:.4f}' for nm, v in zip(VARS, eq)))
    print('All runs start from 1.05 x this state (repository convention).')
    print('=' * 100)

    scenarios = []
    for cs in args.cs:
        sc = run_scenario(cs)
        scenarios.append(sc)
        e = sc['eig']
        eT = f"{e['T_years']:.1f}" if e['T_years'] else 'N/A'
        sC = sc['var_stats']['sC']; sV = sc['var_stats']['sV']
        print(f"\nCS = {cs:+.3f}   (first sV<0 at t_dom = {sc['t_dom']:.1f} y; "
              f"first share>1 at {sc['t_dom_shares']:.1f} y)")
        print(f"  eigenvalue period at reference state: {eT} y (Re = {e['Re']:+.4f})")
        print(f"  sC: {sC['n_cycles']} cycles {sC['cycles']} | T_med {sC['T_med']:.1f} y | "
              f"rise {sC['rise_med']:.1f} y | fall {sC['fall_med']:.1f} y | "
              f"rise share {sC['rise_share_pct']:.0f}%")
        print(f"  sV: {sV['n_cycles']} cycles {sV['cycles']} | T_med {sV['T_med']:.1f} y | "
              f"rise {sV['rise_med']:.1f} y | fall {sV['fall_med']:.1f} y | "
              f"rise share {sV['rise_share_pct']:.0f}%")
        print(f"  r inside t_dom: mean {sc['r_mean']:.4f}, min {sc['r_min']:.4f}, "
              f"end {sc['r_end']:.4f}")

    pulse = None
    ratchet = None
    if not args.no_figs:
        pulse = _run_pulse_episode(eq)
        ratchet = run_ratchet_recoveries()
        _make_figures(scenarios, eq, pulse, ratchet)

    os.makedirs(DATA_DIR, exist_ok=True)
    report = {
        'params': P,
        'note': ('CS-extended model (8.15-CS)-(8.19); CS=0 reproduces the book. '
                 'ge = -sV - CS in (8.15). Metrics are measured inside the economically '
                 'meaningful window t < t_dom = first time sV<0 (repository convention); '
                 't_dom_shares = first time sC>1 or sV>1 (book premise sC+sV<=1) is '
                 'reported separately.'),
        'scenarios': [dict((k, v) for k, v in sc.items() if k not in ('t', 'sol', 'msk'))
                      for sc in scenarios],
    }
    if pulse is not None:
        report['pulse_episode'] = pulse['summary']
    if ratchet is not None:
        report['ratchet_recoveries'] = ratchet['summary']
    def _clean(o):
        """Recursively convert numpy types and NaN to strict-JSON values
        (NaN occurs for variables with no complete cycle in the window,
        e.g. the profit rate r, which only declines within t_dom)."""
        if isinstance(o, dict):
            return {k: _clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_clean(v) for v in o]
        if isinstance(o, (float, np.floating)):
            return None if math.isnan(float(o)) else round(float(o), 6)
        if isinstance(o, (int, np.integer)):
            return int(o)
        return o
    with open(args.out, 'w', encoding='utf-8') as fh:
        json.dump(_clean(report), fh, indent=2, ensure_ascii=False)
    print(f"\nResults written to {args.out}")


def _run_pulse_episode(eq, T_years=150.0):
    """Stylised 'historical episode' experiment: a militant wage-push phase
    followed by a capital offensive, against the constant-CS=0 baseline.

    Profile of CS over model time (interpretation in brackets):
      t in [0, 18):      0.000   - 'normal' struggle level
      t in [18, 30):    +0.020   - a wage explosion at the boom peak
                                   (the late-1960s/early-1970s analogy)
      t in [30, 44):    +0.005   - persistent above-normal push after the
                                   explosion (the high-wage-share 1970s)
      t in [44, 150):   -0.015   - the capital counter-offensive
                                   (the post-1980 analogy)
    Model time is not calibrated to calendar years; the run is illustrative."""
    def cs_profile(t):
        if t < 18:   return 0.0
        if t < 30:   return 0.02
        if t < 44:   return 0.005
        return -0.015
    y0 = 1.05 * eq
    t = np.linspace(0, T_years, int(T_years * 80) + 1)
    sol_p = odeint(lambda y, t: f(y, cs_profile(t)), y0, t, rtol=1e-9, atol=1e-12)
    sol_b = odeint(lambda y, t: f(y, 0.0), y0, t, rtol=1e-9, atol=1e-12)
    prof = np.array([cs_profile(tt) for tt in t])
    # scalar summary (also written to the JSON report)
    def at(sol, tt):
        i = int(tt * 80)
        return float(sol[i, 0])
    m30_55 = (t >= 30) & (t <= 55)
    m55_100 = (t >= 55) & (t <= 100)
    summary = dict(
        r_at=[round(at(s, tt), 4) for s in (sol_p, sol_b) for tt in (18, 30, 44, 60, 80, 100)],
        r_mean_downturn=[float(sol_p[m30_55, 0].mean()), float(sol_b[m30_55, 0].mean())],
        r_min_downturn=[float(sol_p[m30_55, 0].min()), float(sol_b[m30_55, 0].min())],
        r_peak_recovery=[float(sol_p[m55_100, 0].max()), float(sol_b[m55_100, 0].max())],
        sV_min_squeeze=[float(sol_p[m30_55, 1].min()), float(sol_b[m30_55, 1].min())],
        sV_at_60=[float(sol_p[60 * 80, 1]), float(sol_b[60 * 80, 1])],
    )
    print('\nPulse experiment (wage explosion -> offensive vs constant CS=0):')
    print(f"  mean r over the downturn t in [30,55]: pulse {summary['r_mean_downturn'][0]:.4f}"
          f" vs baseline {summary['r_mean_downturn'][1]:.4f}")
    print(f"  min  r over the downturn:              pulse {summary['r_min_downturn'][0]:.4f}"
          f" vs baseline {summary['r_min_downturn'][1]:.4f}")
    print(f"  peak r over the recovery t in [55,100]: pulse {summary['r_peak_recovery'][0]:.4f}"
          f" vs baseline {summary['r_peak_recovery'][1]:.4f}")
    print(f"  min sV over the squeeze t in [30,55]:  pulse {summary['sV_min_squeeze'][0]:.4f}"
          f" vs baseline {summary['sV_min_squeeze'][1]:.4f}")
    return dict(t=t, profile=prof, sol_pulse=sol_p, sol_base=sol_b, summary=summary)


def run_ratchet_recoveries(cs_fixed=0.02, cs_ref=0.0, T_years=220.0):
    """The 'pointless recoveries' run: CS fixed at a positive value for the
    whole integration (a run with CS > 0 keeps producing accumulation cycles -
    every upswing restarts when the hiring pressure sV*r falls back below the
    choke threshold tau, without any restoration of the rate of surplus value),
    while the profit rate ratchets down cycle after cycle.

    Returns the trajectories, the sC events (peaks/troughs) with the profit
    rate at each event, the domain-exit markers and a scalar summary (also
    written to the JSON report). cs_fixed = +0.02 is the 1970s-scale wage push
    held constant; cs_ref = 0 is the book baseline for comparison."""
    eq = equilibrium_closed_form()
    y0 = 1.05 * eq
    t = np.linspace(0, T_years, int(T_years * 80) + 1)
    sol_f = odeint(lambda y, t: f(y, cs_fixed), y0, t, rtol=1e-9, atol=1e-12)
    sol_r = odeint(lambda y, t: f(y, cs_ref), y0, t, rtol=1e-9, atol=1e-12)
    # sC peaks/troughs of the fixed-CS run (skip the initial transient peak)
    pk, tr = peak_trough_indices(t, sol_f[:, 2], 10.0)
    ev = sorted([(float(t[i]), 'P') for i in pk if t[i] > 12]
                + [(float(t[i]), 'T') for i in tr])
    peaks = [e for e in ev if e[1] == 'P']
    troughs = [e for e in ev if e[1] == 'T']
    def r_at(sol, tt):
        i = int(np.searchsorted(t, tt))
        return float(sol[min(i, len(t) - 1), 0])
    neg = np.where(sol_f[:, 1] < 0)[0]
    gt1 = np.where(sol_f[:, 2] > 1.0)[0]
    t_sV0 = float(t[neg[0]]) if len(neg) else None
    t_sC1 = float(t[gt1[0]]) if len(gt1) else None
    summary = dict(
        cs=cs_fixed, cs_ref=cs_ref,
        peak_times=[round(e[0], 1) for e in peaks[:6]],
        r_at_peaks=[round(r_at(sol_f, e[0]), 4) for e in peaks[:6]],
        r_at_peaks_ref=[round(r_at(sol_r, e[0]), 4) for e in peaks[:6]],
        trough_times=[round(e[0], 1) for e in troughs[:6]],
        r_at_troughs=[round(r_at(sol_f, e[0]), 4) for e in troughs[:6]],
        peak_spacings=[round(peaks[i+1][0] - peaks[i][0], 1) for i in range(min(5, len(peaks)-1))],
        t_sC_gt_1=t_sC1, t_sV_lt_0=t_sV0,
        r_120=round(r_at(sol_f, 120.0), 4), r_220=round(r_at(sol_f, 220.0), 4),
        r_120_ref=round(r_at(sol_r, 120.0), 4), r_220_ref=round(r_at(sol_r, 220.0), 4),
    )
    print('\nPointless recoveries at a fixed wage push (CS = %+.3f, book baseline CS = %+.3f):'
          % (cs_fixed, cs_ref))
    print('  sC boom peaks at t = %s with r = %s (CS = 0 at the same dates: %s)'
          % (summary['peak_times'], summary['r_at_peaks'], summary['r_at_peaks_ref']))
    print('  sC troughs  at t = %s with r = %s' % (summary['trough_times'], summary['r_at_troughs']))
    print('  peak-to-peak spacings: %s  |  sC > 1 at t = %s  |  sV < 0 at t = %s'
          % (summary['peak_spacings'], t_sC1, t_sV0))
    print('  r(120) = %.4f (ref %.4f)  r(220) = %.4f (ref %.4f)'
          % (summary['r_120'], summary['r_120_ref'], summary['r_220'], summary['r_220_ref']))
    return dict(t=t, sol_fixed=sol_f, sol_ref=sol_r, peaks=peaks, troughs=troughs,
                summary=summary)


def _make_figures(scenarios, eq, pulse=None, ratchet=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    cs_all = [sc['cs'] for sc in scenarios]
    colors = plt.cm.coolwarm(np.linspace(0.05, 0.95, len(scenarios)))

    # ---------------- Figure 1: trajectories (r, sV, sC) per CS ----------------
    fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
    for sc, col in zip(scenarios, colors):
        msk = sc['msk']
        for ax, idx, nm in zip(axes, [0, 1, 2], ['r', 'sV', 'sC']):
            ax.plot(sc['t'][msk], sc['sol'][msk, idx], lw=1.1, color=col,
                    label=f"CS={sc['cs']:+.2f}")
            ax.set_ylabel(f'{nm}(t)')
            ax.grid(alpha=0.3)
        axes[0].axhline(eq[0], color='k', ls=':', lw=0.8)
    axes[0].legend(fontsize=8, ncol=4, loc='upper right')
    axes[0].set_title('Profit rate r, sV and sC under different class-struggle levels CS\n'
                      '(each trajectory is truncated at its first sV<0, the domain exit)')
    axes[2].set_xlabel('t (years)')
    fig.tight_layout()
    f1 = os.path.join(FIG_DIR, 'class_struggle_trajectories.png')
    fig.savefig(f1, dpi=130); plt.close(fig)

    # ---------------- Figure 2: response of cycle length and slope lengths -------
    sC_stats = [sc['var_stats']['sC'] for sc in scenarios]
    sV_stats = [sc['var_stats']['sV'] for sc in scenarios]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    ax.plot(cs_all, [s['T_med'] for s in sC_stats], 'o-', label='sC (nonlinear)')
    ax.plot(cs_all, [s['T_med'] for s in sV_stats], 's--', label='sV (nonlinear)')
    ax.plot(cs_all, [sc['eig']['T_years'] for sc in scenarios], 'k:',
            label='eigenvalue period at reference state')
    ax.set_xlabel('CS'); ax.set_ylabel('median cycle length (years)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title('Cycle length vs class struggle')
    ax = axes[0, 1]
    ax.plot(cs_all, [s['rise_med'] for s in sC_stats], 'o-', label='sC rise')
    ax.plot(cs_all, [s['fall_med'] for s in sC_stats], 'o--', label='sC fall')
    ax.plot(cs_all, [s['rise_med'] for s in sV_stats], 's-', label='sV rise')
    ax.plot(cs_all, [s['fall_med'] for s in sV_stats], 's--', label='sV fall')
    ax.set_xlabel('CS'); ax.set_ylabel('median duration (years)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title('Upswing (rise) and downswing (fall) lengths')
    ax = axes[1, 0]
    ax.plot(cs_all, [s['rise_share_pct'] for s in sC_stats], 'o-', label='sC')
    ax.plot(cs_all, [s['rise_share_pct'] for s in sV_stats], 's--', label='sV')
    ax.axhline(50, color='k', ls=':', lw=0.8)
    ax.set_xlabel('CS'); ax.set_ylabel('rise share of the cycle (%)')
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax.set_title('Asymmetry: share of the cycle spent in the upswing')
    ax = axes[1, 1]
    ax.plot(cs_all, [sc['r_mean'] for sc in scenarios], 'o-', label='mean r (in-domain)')
    ax.plot(cs_all, [sc['r_end'] for sc in scenarios], 's--', label='r at t_dom')
    ax.set_xlabel('CS')
    ax.set_ylabel('profit rate')
    ax.legend(fontsize=8, loc='center right'); ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(cs_all, [sc['t_dom'] for sc in scenarios], '^:', color='green')
    ax2.set_ylabel('t_dom (years, green)')
    ax.set_title('Profitability and the duration of the meaningful domain')
    fig.suptitle('Effect of the class-struggle level CS on the long-wave cycle '
                 '(book parameters, n = 0.015)', y=0.99)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    f2 = os.path.join(FIG_DIR, 'class_struggle_cycle_metrics.png')
    fig.savefig(f2, dpi=130); plt.close(fig)

    print(f"\nFigures saved: {f1}\n                {f2}")

    # ---------------- Figure 3: the pulse / epoch-switch experiment -----------
    if pulse is not None:
        t = pulse['t']
        fig, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        # shade the three CS regimes
        regimes = [(18, 30, '#d62728', 'wage explosion, CS=+0.02'),
                   (30, 44, '#ff9896', 'persistent push, CS=+0.005'),
                   (44, 150, '#1f77b4', 'capital offensive, CS=-0.015')]
        for a, b, col, lbl in regimes:
            for ax in axes:
                ax.axvspan(a, b, color=col, alpha=0.10)
            axes[0].text((a+b)/2, 1.01, lbl, ha='center', va='bottom', fontsize=7.5,
                         transform=axes[0].get_xaxis_transform(), color=col)
        # deviations of the pulse run from the CS=0 baseline (the effects are
        # small - the model's distribution channel carries the weight a2=0.01 -
        # so deviations, not levels, are plotted)
        for idx, ylab, unit in [(0, 'deviation of r from baseline (%)', '%'),
                                (1, 'deviation of sV from baseline', ''),
                                (2, 'deviation of sC from baseline', '')]:
            ax = axes[idx]
            d = pulse['sol_pulse'][:, idx] - pulse['sol_base'][:, idx]
            if idx == 0:
                d = 100 * d / pulse['sol_base'][:, idx]
            ax.plot(t, d, lw=1.2, color='#1f77b4')
            ax.axhline(0, color='#555555', lw=0.8)
            ax.set_ylabel(ylab); ax.grid(alpha=0.3)
            imax = int(np.argmax(np.abs(d)))
            ax.annotate(f'max |dev| {d[imax]:+.2f}{unit} at t={t[imax]:.0f}',
                        xy=(t[imax], d[imax]), xytext=(t[imax]-55, d[imax] + 0.25*np.std(d) + 1e-9),
                        fontsize=7.5, arrowprops=dict(arrowstyle='->', lw=0.7))
        axes[0].set_title('Stylised historical episode (pulse run minus CS = 0 baseline):\n'
                          'militant wage-push phase at t = 18-44, capital offensive t > 44. '
                          'The class-struggle episode moves the profit rate by at most ~1%\n'
                          'because in the book calibration the wage-push channel enters with the '
                          'small weight a2 = 0.01 (model time is not calibrated to calendar years)')
        axes[2].set_xlabel('t (years)')
        fig.tight_layout()
        f3 = os.path.join(FIG_DIR, 'class_struggle_pulse_episode.png')
        fig.savefig(f3, dpi=130); plt.close(fig)
        print(f"                {f3}")

    # -------- Figure 4: pointless recoveries at a fixed CS > 0 (the ratchet) ----
    if ratchet is not None:
        t = ratchet['t']
        s = ratchet['summary']
        sol = ratchet['sol_fixed']
        hr = sol[:, 1] * sol[:, 0]          # hiring pressure sV*r
        tau = sol[:, 4]
        rec_mask = hr < tau                 # sC rising: the "recovery" phases
        fig, axes = plt.subplots(4, 1, figsize=(11, 13.5), sharex=True)
        # shade recoveries (hr < tau) in light green on every panel
        for ax in axes:
            ax.fill_between(t, 0, 1, where=rec_mask, transform=ax.get_xaxis_transform(),
                            color='#2ca02c', alpha=0.07)
        # -- panel 1: the ratchet of the profit rate
        ax = axes[0]
        ax.plot(t, ratchet['sol_ref'][:, 0], lw=1.1, color='#888888', ls='--',
                label=f'baseline CS = {s["cs_ref"]:+.2f} (book)')
        ax.plot(t, sol[:, 0], lw=1.5, color='#d62728',
                label=f'fixed CS = {s["cs"]:+.2f} (the push never ends)')
        for i, (tt, rv) in enumerate(zip(s['peak_times'], s['r_at_peaks'])):
            ax.plot(tt, rv, 'o', ms=5, color='#d62728')
            off = 0.006 if i % 2 == 0 else -0.011
            ax.annotate(f'r = {rv:.3f}', xy=(tt, rv), xytext=(tt - 6, rv + off),
                        fontsize=8, color='#d62728')
        ax.plot(s['trough_times'][0], s['r_at_troughs'][0], 's', ms=5, color='#d62728')
        ax.set_ylabel('r(t)'); ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.text(0.02, 0.97, 'each boom peak is lower than the last: the recoveries restart,\n'
                'but every restart begins from a lower profit rate (the ratchet)',
                transform=ax.transAxes, va='top', fontsize=8.5,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85, ec='#888888', lw=0.5))
        # -- panel 2: the accumulation wave keeps cycling
        ax = axes[1]
        ax.plot(t, sol[:, 2], lw=1.3, color='#1f77b4')
        for i, e in enumerate(ratchet['peaks'][:5]):
            ax.annotate(f'P{i+1}', xy=(e[0], sol[int(np.searchsorted(t, e[0])), 2]),
                        xytext=(e[0] - 3, 0), textcoords='offset points', fontsize=8)
        for i, e in enumerate(ratchet['troughs'][:5]):
            ax.annotate(f'T{i+1}', xy=(e[0], sol[int(np.searchsorted(t, e[0])), 2]),
                        xytext=(e[0] - 3, -12), textcoords='offset points', fontsize=8)
        ax.axvline(s['t_sC_gt_1'], color='#555555', ls=':', lw=1.2)
        ax.text(s['t_sC_gt_1'] + 1, ax.get_ylim()[1]*0.97, 'sC > 1\n(share above the whole\n'
                'surplus value: domain ends)', fontsize=7.5, va='top')
        ax.set_ylabel('sC(t)'); ax.grid(alpha=0.3)
        ax.text(0.02, 0.97, 'the mechanisation wave continues at fixed CS > 0: peaks P1-P4 at '
                f't = {s["peak_times"][0]:.0f}, {s["peak_times"][1]:.0f}, '
                f'{s["peak_times"][2]:.0f} and {s["peak_times"][3]:.0f} y\n'
                '(spacings ~47-63 y, close to the CS = 0 cadence): no restoration of the '
                'surplus-value rate is needed for the boom to restart',
                transform=ax.transAxes, va='top', fontsize=8.5,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85, ec='#888888', lw=0.5))
        # -- panel 3: the restart mechanism (hiring pressure vs the choke threshold)
        ax = axes[2]
        ax.plot(t, hr, lw=1.2, color='#9467bd', label='hiring pressure  h = sV*r')
        ax.plot(t, tau, lw=1.2, color='#2ca02c', label='choke threshold  tau(t)')
        ax.set_ylabel('h = sV*r, tau'); ax.grid(alpha=0.3)
        ax.legend(fontsize=8, loc='upper right')
        ax.text(0.02, 0.97, 'green shading = h < tau, i.e. sC rising (the "recovery" phases).\n'
                'Each recovery starts when h falls back below tau - it does NOT wait for the\n'
                'rate of surplus value to be restored (g_e < 0 all along while h > -CS)',
                transform=ax.transAxes, va='top', fontsize=8.5,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85, ec='#888888', lw=0.5))
        # -- panel 4: the hiring share and how the run ends
        ax = axes[3]
        ax.plot(t, sol[:, 1], lw=1.2, color='#ff7f0e')
        ax.axhline(0, color='k', lw=0.7)
        ax.axvline(s['t_sV_lt_0'], color='#d62728', ls='--', lw=1.2)
        ax.text(s['t_sV_lt_0'] + 1, ax.get_ylim()[0], 'sV < 0 at t = %.0f y:\nthe purge arrives\n'
                'only late (CS below the\npurge threshold)' % s['t_sV_lt_0'], fontsize=7.5, va='bottom')
        ax.set_ylabel('sV(t)'); ax.grid(alpha=0.3)
        ax.text(0.02, 0.97, 'the hiring share stays positive for ~133 y: the reserve-army purge is '
                'delayed.\nr(120) = %.3f, r(220) = %.3f (baseline %.3f / %.3f): the profit slide '
                'continues -\nthe profitless succession has no internal end; only the regime change '
                '(CS < 0) restores profit'
                % (s['r_120'], s['r_220'], s['r_120_ref'], s['r_220_ref']),
                transform=ax.transAxes, va='top', fontsize=8.5,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85, ec='#888888', lw=0.5))
        axes[3].set_xlabel('t (years)')
        fig.suptitle('Pointless recoveries at a fixed wage push: CS = %+.2f held constant over the '
                     'whole run (book parameters, n = 0.015)\n'
                     'the accumulation wave keeps cycling while the profit rate ratchets down - '
                     'model time, not calibrated to calendar years' % s['cs'], y=0.995)
        fig.tight_layout(rect=[0, 0, 1, 0.975])
        f4 = os.path.join(FIG_DIR, 'class_struggle_ratchet_recoveries.png')
        fig.savefig(f4, dpi=130); plt.close(fig)
        print(f"                {f4}")


if __name__ == '__main__':
    main()
