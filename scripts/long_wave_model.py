# -*- coding: utf-8 -*-
"""
Author: sacanlee

Reproduction of the long-wave (Kondratiev) model and scenario analysis over the
working-age population growth rate n.

Reference: Chatzarakis, Tsaliki & Tsoulfidis (2022), "Economic Growth and Long
Cycles: A Classical Political Economy Approach", Routledge, Ch. 8, Section
8.3.3 "Complete Model".

Model (Eqs. 8.15-8.19 of the book):
    r'    = -a1*(sC - sV)*r**2 + a2*(delta + tau - sV)*r      (rate of profit)
    sV'   = sC*r - delta - tau - n                             (variable-capital share)
    sC'   = sC*(tau - sV*r)                                    (constant-capital share)
    delta'= b0 + delta*(gw - sC*r)                             (devaluation rate)
    tau'  = b1*Erf(delta - b2)                                 (technical change rate)

Equilibrium (printed closed form, Eq. 8.25):
    r*     = a2*(b0 - b2*(b2 - gw + n)) / (a1*(b2*(b2+n) - a2*(b0 + b2*(gw - n))))
    sV*    = (a1/a2)*(b2 + n) - b0/b2 - gw + n
    sC*    = (b0 + b2*gw)*(a1*(b2*(b2+n) - a2*(b0 + b2*(gw - n)))) / (a2*b2*(b0 - b2*(b2 - gw + n)))
    delta* = b2
    tau*   = b0/b2 - b2 + gw - n

Book parameters: a1=0.04, a2=0.01, b0=0.004, b1=0.005, b2=0.05, n=0.015, gw=0.03.
The book claims the equilibrium {r*, sV*, sC*, delta*, tau*} = {0.0947, 0.1188,
0.8708, 0.05, 0.045} and a cycle of "about 50 years"; its own printed matrix
eigenvalues imply 2*pi/0.0947 = 66.3 years (see script output and docs).

Usage:
    python long_wave_model.py                 # all self-checks + default 4 n scenarios
    python long_wave_model.py --ns 0.015 0 -0.008 -0.015   # custom scenarios
    python long_wave_model.py --ns 0.02 0.015 0.005        # positive n range (nonlinear reference)
    python long_wave_model.py --a1 0.04 --n -0.010         # custom model parameters
    python long_wave_model.py --cycles-fig figures/ode_solution_book_params.png
                                                 # figure: nonlinear solution cycles
                                                 # at the given (book) parameters
"""
import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks
import sympy as sp
import math, json, argparse, os, sys

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')
    sys.stderr.reconfigure(encoding='utf-8')

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..'))
DATA_DIR = os.path.join(REPO_ROOT, 'data')
FIG_DIR = os.path.join(REPO_ROOT, 'figures')

# ----------------------------------------------------------------------
# Symbolic setup: analytical Jacobian (equations as printed on paper)
# ----------------------------------------------------------------------
_syms_r, _syms_sV, _syms_sC, _syms_d, _syms_t = sp.symbols('r sV sC d tau')
_a1, _a2, _b0, _b1, _b2, _gw, _n = sp.symbols('a1 a2 b0 b1 b2 gw n')

def _sym_jacobian():
    """Analytical Jacobian d(f)/d(x) at (r,sV,sC,d,tau). Erf is linearised at
    the equilibrium delta=b2, with erf'(0) = 2/sqrt(pi)."""
    r, sV, sC, d, t = _syms_r, _syms_sV, _syms_sC, _syms_d, _syms_t
    rdot = -_a1*(sC - sV)*r**2 + _a2*(d + t - sV)*r
    sVdot = sC*r - d - t - _n
    sCdot = sC*(t - sV*r)
    ddot = _b0 + d*(_gw - sC*r)
    tdot = _b1*(2/sp.sqrt(sp.pi))*(d - _b2)
    J = sp.Matrix([[sp.diff(f, v) for v in (r, sV, sC, d, t)]
                   for f in (rdot, sVdot, sCdot, ddot, tdot)])
    return sp.lambdify((r, sV, sC, d, t, _a1, _a2, _b0, _b1, _b2, _gw, _n), J, 'numpy')

JAC_F = _sym_jacobian()

# ----------------------------------------------------------------------
# Model core
# ----------------------------------------------------------------------
class Model:
    def __init__(self, a1=0.04, a2=0.01, b0=0.004, b1=0.005, b2=0.05, gw=0.03, n=0.015,
                 erf_deriv=True):
        self.p = dict(a1=a1, a2=a2, b0=b0, b1=b1, b2=b2, gw=gw, n=n)
        self.erf_deriv = erf_deriv

    def f(self, y):
        r, sV, sC, d, tau = y
        a1, a2, b0, b1, b2, gw, n = [self.p[k] for k in ('a1','a2','b0','b1','b2','gw','n')]
        return np.array([-a1*(sC - sV)*r**2 + a2*(d + tau - sV)*r,
                         sC*r - d - tau - n,
                         sC*(tau - sV*r),
                         b0 + d*(gw - sC*r),
                         b1*math.erf(d - b2)])

    # ---------------- Equilibrium ----------------
    def equilibrium_closed_form(self):
        """Closed-form solution of the book's Eq. (8.25). Note: the equilibrium
        values printed in the book {0.0947,0.1188,0.8708,0.05,0.045} do NOT match
        this formula evaluated at the book parameters (see book_selfcheck)."""
        a1, a2, b0, b2, gw, n = [self.p[k] for k in ('a1','a2','b0','b2','gw','n')]
        D = b0 - b2*(b2 - gw + n)
        N = b2*(b2 + n) - a2*(b0 + b2*(gw - n))
        r = a2*D/(a1*N)
        sV = (a1/a2)*(b2 + n) - b0/b2 - gw + n
        sC = (b0 + b2*gw)*a1*N/(a2*b2*D)
        d = b2
        tau = b0/b2 - b2 + gw - n
        return np.array([r, sV, sC, d, tau])

    def equilibrium_internal(self):
        """Equilibrium derived from the equations themselves (numerical):
        (8.19)->delta=b2, (8.18)->sC*r=gw+b0/b2, (8.16)->tau, (8.17)->sV=tau/r,
        and (8.15) then yields r* = a2*tau/(a2*(b2+tau)-a1*(b2+n))."""
        a1, a2, b0, b2, gw, n = [self.p[k] for k in ('a1','a2','b0','b2','gw','n')]
        d = b2
        tau = b0/b2 - b2 + gw - n
        r = a2*tau/(a2*(b2 + tau) - a1*(b2 + n))
        sV = tau/r
        sC = (gw + b0/b2)/r
        return np.array([r, sV, sC, d, tau])

    def jacobian(self, x):
        e = (x[0], x[1], x[2], x[3], x[4])
        return JAC_F(*e, *[self.p[k] for k in ('a1','a2','b0','b1','b2','gw','n')])

    # ---------------- Eigen analysis ----------------
    def eig_analysis(self, x):
        J = self.jacobian(x)
        w = np.linalg.eigvals(J)
        return J, w

    @staticmethod
    def describe_eigs(w):
        out = []
        for lam in sorted(w, key=lambda z: -z.imag):
            im = abs(lam.imag)
            T = 2*np.pi/im if im > 1e-10 else float('inf')
            out.append(dict(lam=f"{lam.real:+.5f}{lam.imag:+.5f}i",
                            T_years=(round(T, 1) if T < 1e6 else None)))
        return out

    # ---------------- Simulation ----------------
    def simulate(self, y0, T_years=400, npts=40001):
        t = np.linspace(0, T_years, npts)
        sol = odeint(lambda y, t: self.f(y), y0, t, rtol=1e-9, atol=1e-12)
        return t, sol

    def linear_sim(self, x0, J, T_years=400, npts=40001):
        """Exact solution of the linearised system x' = J x (matrix exponential)."""
        from scipy.linalg import expm
        t = np.linspace(0, T_years, npts)
        sol = np.array([expm(J*ti) @ x0 for ti in t])
        return t, sol

# ----------------------------------------------------------------------
# Book self-check (verify the script against the numbers printed in the book)
# ----------------------------------------------------------------------
BOOK_EQUILIBRIUM = [0.0947, 0.1188, 0.8708, 0.05, 0.045]     # book p. 232
BOOK_MATRIX = np.array([
    [0.0,      -0.0825,  -0.1034,  0.2177,   0.0],
    [0.1263,    0.0,      1.1611,  -1.0,     0.0],
    [-0.00027,  0.00122, -0.0021,  0.00095,  0.00095],
    [0.0,       0.0,      0.0,      0.0,     0.00226],
    [-0.0063,   0.0,     -0.05806, 0.0,     -0.08],
])
BOOK_EIGS = [ -0.07904, -0.0002+0.0947j, -0.0002-0.0947j, -0.0014+0.0001j, -0.0014-0.0001j ]

def book_selfcheck():
    """Consistency checks of the printed values, as a table."""
    m = Model()
    rows = []

    # 1. Claimed equilibrium vs (8.25) closed form
    cf = m.equilibrium_closed_form()
    names = ['r*', 'sV*', 'sC*', 'delta*', 'tau*']
    for k, (nm, v_book, v_cf) in enumerate(zip(names, BOOK_EQUILIBRIUM, cf)):
        rows.append(dict(check=f"Equilibrium {nm}",
                         book=v_book, formula=round(float(v_cf), 5),
                         match=abs(v_book - v_cf) < 0.0005))

    # 2. Eigenvalues of the printed matrix vs eigenvalues printed in the book
    w = np.linalg.eigvals(BOOK_MATRIX)
    w_sorted = sorted(w, key=lambda z: (-z.imag))
    eigs_book = sorted(BOOK_EIGS, key=lambda z: (-z.imag))
    lam_pairs = [(round(eigs_book[i].real, 5), round(w_sorted[i].real, 5),
                  round(eigs_book[i].imag, 5), round(w_sorted[i].imag, 5))
                 for i in range(5)]
    T_book_eig = 2*np.pi/0.0947
    rows.append(dict(check="Dominant complex pair (matrix -> eigenvalues)",
                     book="−0.0002±0.0947i", formula="−0.0002±0.0947i",
                     match=all(abs(a-b) < 2e-3 for a, b in
                               [(lam_pairs[0][0], lam_pairs[0][1]),
                                (lam_pairs[0][2], lam_pairs[0][3])])))
    rows.append(dict(check="Dominant period from the book's linearisation",
                     book='"about 50 years" (text/Fig. 8.13)',
                     formula=round(T_book_eig, 1),
                     match=False))

    # 3. Residuals of the claimed equilibrium in the printed equations
    res = m.f(BOOK_EQUILIBRIUM)
    rows.append(dict(check="Residuals at the claimed equilibrium [r', sV', sC', delta', tau']",
                     book="-", formula="[" + ", ".join(f"{v:+.5f}" for v in res) + "]",
                     match=max(abs(res)) < 1e-3))
    return rows

def book_matrix_eigsonly():
    w = np.linalg.eigvals(BOOK_MATRIX)
    return sorted(w, key=lambda z: -z.imag)

def _peak_spacings(t, yv):
    """Times and lengths of the completed peak-to-peak intervals of a series.
    Returns (spacing_start_times, spacings_in_years)."""
    from scipy.signal import find_peaks
    pk, _ = find_peaks(yv, distance=max(2, len(t) // 200))
    if len(pk) < 2:
        return [], []
    t_pk = t[pk]
    return t_pk[1:], list(np.diff(t_pk))

def plot_solution_cycles(m, out_png, T_years=300.0, T_long=600.0):
    """Figure: the cycles of the nonlinear ODE solution at the given parameters.

    Upper block: r, sV, sC, delta, tau over T_years (default 300) years, with
    the measured peak-to-peak spacings annotated per panel. Initial values are
    the closed-form equilibrium (8.25) times 1.05 (same convention as in
    run_scenario).

    Bottom strip: peak-to-peak spacing of every completed cycle vs time from a
    T_long (600 y) integration. The point of the figure is that the spacing
    does NOT settle at a stable value: amplitudes grow and the period keeps
    lengthening (sC 43 -> 83 y over 300 y; sV/delta/tau up to ~140 y by 600 y)
    because the closed-form equilibrium is an unstable focus (Re = +0.0212).
    The trajectory never reaches a limit cycle, so "about 50 years" is only a
    median over the first cycles - see docs/model_validation_n_scenarios.md
    and docs/Technical_Report_Ch8_Model_Errata.docx (Finding 6)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    names = ['r', 'sV', 'sC', 'delta', 'tau']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#d62728']
    eq = m.equilibrium_closed_form()
    twin = T_years                      # variable-panel window (>= 300 y)
    t, sol = m.simulate(eq * 1.05, T_years=twin)
    tL, solL = m.simulate(eq * 1.05, T_years=T_long)   # long run for the spacing strip

    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(10, 15.6))
    gs = gridspec.GridSpec(6, 1, height_ratios=[2.2, 2.2, 2.2, 2.2, 2.2, 1.5],
                           hspace=0.38, top=0.94, bottom=0.075, left=0.09, right=0.97)
    axes = [fig.add_subplot(g) for g in gs]

    for ax, idx, nm, col in zip(axes[:5], range(5), names, colors):
        ax.plot(t, sol[:, idx], lw=1.0, color=col)
        ax.axhline(eq[idx], color='#666666', ls='--', lw=0.9,
                   label=f"equilibrium (8.25): {nm}* = {eq[idx]:.4f}")
        starts, sp = _peak_spacings(t, sol[:, idx])
        if sp:
            med = float(np.median(sp))
            txt = (f"{len(sp)} complete cycles in {twin:.0f} y  |  peak-to-peak spacing "
                   f"{sp[0]:.0f} -> {sp[-1]:.0f} y  |  median {med:.0f} y\n"
                   f"-> spacing lengthens, no stable period")
        else:
            txt = (f"no complete cycle in {twin:.0f} y: r decays from {sol[0, idx]:.2f} "
                   f"towards {sol[-1, idx]:.3f}\n-> tendency of the rate of profit to fall")
        ax.text(0.015, 0.94, txt, transform=ax.transAxes, fontsize=8.6, va='top',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.88,
                          ec='#888888', lw=0.5))
        ax.set_ylabel(f'{nm}(t)', fontsize=10)
        ax.legend(loc='lower right', fontsize=8, framealpha=0.9)
        ax.grid(alpha=0.3)
        if idx < 4:
            ax.set_xticklabels([])
    axes[4].set_xlabel('t (years)', fontsize=10)

    # Bottom strip: spacing of every completed cycle vs time (600-y integration).
    axs = axes[5]
    for idx, nm, col in zip([1, 2, 3, 4], ['sV', 'sC', 'delta', 'tau'], colors[1:]):
        starts, sp = _peak_spacings(tL, solL[:, idx])
        # After ~350 y sC stops oscillating and explodes monotonically (the
        # spacing ~757 y interval is an artefact of that regime); keep sC only
        # while it still completes cycles within the 600-y window.
        if nm == 'sC':
            keep = [s for s in sp if s < 200]
            starts, sp = starts[:len(keep)], keep
        if sp:
            axs.plot(starts, sp, 'o-', ms=4, lw=1.0, color=col, label=nm)
    axs.axhline(50, color='#555555', ls='--', lw=1.0)
    axs.text(1.0, 50, '  book: "about 50 years"', va='center', ha='right',
             fontsize=8, color='#555555', transform=axs.get_yaxis_transform())
    axs.set_xlim(0, T_long)
    axs.set_ylim(0, 170)
    axs.set_xlabel('time at the end of each cycle (years)', fontsize=9)
    axs.set_ylabel('peak-to-peak\nspacing (y)', fontsize=9)
    axs.legend(loc='upper left', fontsize=8, ncol=4, framealpha=0.9)
    axs.grid(alpha=0.3)
    axs.set_title('Completed-cycle spacing over a 600-year integration: '
                  'the period keeps lengthening and never stabilises '
                  '(sC stops cycling after ~330 y and explodes monotonically)',
                  fontsize=9.5, loc='left')

    p = m.p
    fig.suptitle("Nonlinear solution of the long-wave model (8.15)-(8.19) at the book parameters\n"
                 f"a1={p['a1']}, a2={p['a2']}, b0={p['b0']}, b1={p['b1']}, b2={p['b2']}, "
                 f"gw={p['gw']}, n={p['n']}   |   initial values = closed-form equilibrium "
                 f"(8.25) x 1.05, window {twin:.0f} y",
                 fontsize=11)
    fig.text(0.02, 0.014,
             'Why no stable period: the closed-form equilibrium (8.25) is an unstable focus '
             '(dominant eigenvalue Re = +0.0212), so the orbit spirals out - amplitudes grow\n'
             'about x3 per cycle and the peak spacing lengthens (sC: 43 -> 83 y over 300 y; '
             'up to ~140 y by 600 y).  After ~300 y sV turns negative and the model leaves its\n'
             'economic domain, so longer integrations do not reveal a limit cycle either.  '
             'Eigenvalue periods: 41.7 y at (8.25); 66.3 y (2*pi/0.0947) for the book\'s printed '
             'matrix -\nsee docs/model_validation_n_scenarios.md.',
             ha='left', va='bottom', fontsize=8.2, color='#333333')
    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)

# ----------------------------------------------------------------------
# Rise/fall segment measurement
# ----------------------------------------------------------------------
def cycle_segments(t, y, pick_index=0):
    """Detect peaks/troughs in a time series; return
    (period_med, rise_med, fall_med, n_cycles, win_med)
    - period_med: median peak-to-peak period (true cycle length)
    - rise_med / fall_med: median trough->peak (rise) / peak->trough (fall) duration
    - n_cycles: number of complete cycles; win_med: median adjacent-event spacing (~T/2)"""
    yv = y[:, pick_index]
    pk, _ = find_peaks(yv, distance=len(t)//80)
    tr, _ = find_peaks(-yv, distance=len(t)//80)
    events = sorted([(int(i), 'pk') for i in pk] + [(int(i), 'tr') for i in tr])
    rises, falls, win = [], [], []
    for i in range(1, len(events)):
        dt = t[events[i][0]] - t[events[i-1][0]]
        win.append(dt)
        if events[i][1] == 'pk' and events[i-1][1] == 'tr':
            rises.append(dt)
        elif events[i][1] == 'tr' and events[i-1][1] == 'pk':
            falls.append(dt)
    periods = [ri + fa for ri, fa in zip(rises, falls)]
    def med(x): return float(np.median(x)) if len(x) else float('nan')
    return (med(periods) if periods else float('nan'),
            med(rises), med(falls), len(periods),
            med(win) if win else float('nan'))

# ----------------------------------------------------------------------
# Main flow
# ----------------------------------------------------------------------
def run_scenario(m, n, linear=True, T_years=300):
    """Return equilibrium, Jacobian, eigenvalues, period, and rise/fall segments
    for a given n."""
    m.p['n'] = n
    eq = m.equilibrium_closed_form()
    J, w = m.eig_analysis(eq)
    desc = m.describe_eigs(w)
    dom = None
    for d0 in desc:
        if d0['T_years'] is not None:
            dom = d0
            break
    # Initial perturbation: ±5% of the equilibrium on all five variables
    # (transparent, reproducible), so that all modes (incl. fast-decaying real
    # modes) are excited - the rise/fall asymmetry comes from mode mixing.
    x0 = 0.05*eq
    t, sol = m.linear_sim(x0, J, T_years=T_years)
    segs = cycle_segments(t, sol, pick_index=0)
    # Nonlinear system: window of 6 dominant periods (for positive n the
    # closed-form equilibrium drifts, so a longer window is needed).
    T0 = dom['T_years'] if dom else 50.0
    nonlin_t, nonlin_sol = m.simulate(eq + x0, T_years=min(T_years, 6*T0))
    segs_nl_r = cycle_segments(nonlin_t, nonlin_sol, pick_index=0)     # r(t)
    segs_nl_sC = cycle_segments(nonlin_t, nonlin_sol, pick_index=2)    # sC(t)
    return dict(n=n, eq=[round(float(v), 5) for v in eq], J=J, eigs=desc,
                dominant=dom, segments=segs, segments_nonlin=segs_nl_r,
                segments_nonlin_sC=segs_nl_sC,
                sim_t=t, sim=sol, nonlin_t=nonlin_t, nonlin=nonlin_sol)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--a1', type=float, default=0.04)
    ap.add_argument('--a2', type=float, default=0.01)
    ap.add_argument('--b0', type=float, default=0.004)
    ap.add_argument('--b1', type=float, default=0.005)
    ap.add_argument('--b2', type=float, default=0.05)
    ap.add_argument('--gw', type=float, default=0.03)
    ap.add_argument('--n', type=float, default=0.015)
    ap.add_argument('--ns', type=float, nargs='+', default=[0.015, 0.0, -0.008, -0.015],
                    help='list of n scenarios (default: book baseline + 0 + -0.008 + -0.015)')
    ap.add_argument('--out', default=os.path.join(DATA_DIR, 'results.json'))
    ap.add_argument('--fig', default=os.path.join(FIG_DIR, 'cycles_by_n.png'))
    ap.add_argument('--cycles-fig', default='',
                    help='optional: PNG for the nonlinear-solution cycles figure at the '
                         'given model parameters (e.g. figures/ode_solution_book_params.png)')
    args = ap.parse_args()

    m = Model(a1=args.a1, a2=args.a2, b0=args.b0, b1=args.b1, b2=args.b2, gw=args.gw,
              n=args.n)

    report = {'params': {k: v for k, v in m.p.items()}, 'selfcheck': book_selfcheck(),
              'scenarios': []}

    print("=" * 78)
    print("A. Book self-check (book values | recomputed values | match)")
    print("=" * 78)
    for r0 in book_selfcheck():
        print(f"{r0['check']:<46s}  {str(r0['book']):<22s}  {str(r0['formula']):<24s}  {r0['match']}")
    wm = book_matrix_eigsonly()
    print(f"\nEigenvalues of the printed matrix (recomputed): " +
          ", ".join(f"{z.real:+.5f}{z.imag:+.5f}i" for z in wm)
          + f"  -> dominant period {2*np.pi/wm[0].imag:.1f} years")

    print("\n" + "=" * 78)
    print("B. n scenarios: equilibrium (8.25 closed form) + analytical Jacobian"
          " eigenvalues / period / rise and fall")
    print("=" * 78)
    scenarios = []
    for n in args.ns:
        sc = run_scenario(m, n)
        scenarios.append(sc)
        report['scenarios'].append({k: v for k, v in sc.items() if k not in ('J','sim_t','sim','nonlin_t','nonlin')})
        e = sc['eq']; d0 = sc['dominant']
        Tstr = f"{d0['T_years']:.1f}" if d0 else "N/A"
        seg = sc['segments']
        seg_nl = sc['segments_nonlin']
        print(f"\nn = {n:+.3f}")
        print(f"  Equilibrium: r*={e[0]:.4f}  sV*={e[1]:.4f}  sC*={e[2]:.4f}  delta*={e[3]:.4f}  tau*={e[4]:.4f}")
        dmod = d0['lam'] if d0 else 'N/A'
        print(f"  Dominant complex eigenvalue: {dmod}  -> eigenvalue period T = {Tstr} years")
        if seg[0] == seg[0]:
            print(f"  Linear system [r(t) peaks/troughs, {seg[3]} cycles]: peak-to-peak {seg[0]:.1f} y | "
                  f"rise {seg[1]:.1f} y | fall {seg[2]:.1f} y | rise share {100*seg[1]/(seg[1]+seg[2]):.0f}%")
        if seg_nl[0] == seg_nl[0]:
            print(f"  Nonlinear system [r(t) peaks/troughs, {seg_nl[3]} cycles]: peak-to-peak {seg_nl[0]:.1f} y | "
                  f"rise {seg_nl[1]:.1f} y ({100*seg_nl[1]/(seg_nl[1]+seg_nl[2]):.0f}%) | "
                  f"fall {seg_nl[2]:.1f} y ({100*seg_nl[2]/(seg_nl[1]+seg_nl[2]):.0f}%)")
        seg_nl2 = sc['segments_nonlin_sC']
        if seg_nl2[0] == seg_nl2[0]:
            print(f"  Nonlinear system [sC(t) peaks/troughs, {seg_nl2[3]} cycles]: peak-to-peak {seg_nl2[0]:.1f} y | "
                  f"rise {seg_nl2[1]:.1f} y ({100*seg_nl2[1]/(seg_nl2[1]+seg_nl2[2]):.0f}%) | "
                  f"fall {seg_nl2[2]:.1f} y ({100*seg_nl2[2]/(seg_nl2[1]+seg_nl2[2]):.0f}%)")

    # Figure: deviation of the profit rate r(t) in the linearised system per scenario
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(len(scenarios), 1, figsize=(9, 2.2*len(scenarios)), sharex=False)
        if len(scenarios) == 1:
            axes = [axes]
        for ax, sc in zip(axes, scenarios):
            T0 = sc['dominant']['T_years']
            win = 2.5*T0
            tt = sc['sim_t']; yy = sc['sim'][:, 0]
            mask = tt <= win if not win else tt <= win
            ax.plot(tt[mask], yy[mask], lw=0.9,
                    label=f"n = {sc['n']:+.3f}   T={T0:.1f}y  "
                          f"rise={sc['segments'][1]:.1f}y / fall={sc['segments'][2]:.1f}y")
            ax.axhline(0, color='k', lw=0.4)
            ax.set_ylabel('Δr(t)')
            ax.legend(loc='upper right', fontsize=8)
            ax.grid(alpha=0.3)
        axes[-1].set_xlabel('t (years)')
        fig.suptitle('Long-wave model (8.15-8.19): linearized deviation of profit rate r(t), '
                     'perturbation = +5%*equilibrium', fontsize=9)
        fig.tight_layout()
        os.makedirs(FIG_DIR, exist_ok=True)
        fig.savefig(args.fig, dpi=130)
        print(f"\nFigure saved: {args.fig}")
    except Exception as e:
        print('Plotting skipped:', e)

    if args.cycles_fig:
        plot_solution_cycles(m, args.cycles_fig)
        print(f"Figure saved: {args.cycles_fig}")

    os.makedirs(DATA_DIR, exist_ok=True)
    with open(args.out, 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else str(o))
    print(f"\nResults written to {args.out}")

if __name__ == '__main__':
    main()
