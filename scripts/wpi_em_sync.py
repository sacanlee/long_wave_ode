# -*- coding: utf-8 -*-
"""
Cycle-synchronisation analysis of the emerging-economy gold-denominated WPI
(PPP-corrected, see wpi_em_build.py) against the United States, plus the
combined figure.

Data:
    data/wpi_em_annual.json              EM goldden series (2010 = 100)
    data/wpi_annual_1800_2026.json       US gold-denominated WPI (rebased 2010=100)

Method ("sync start year"):
    - both series are logged and smoothed with the centred 9-yr moving
      average (as in the WPI cycle work elsewhere in this repository);
    - rolling Pearson correlation over a 20-year window ending at year t
      (>=16 valid years required);
    - the country is "synchronised from year T" if corr(t) >= 0.5 for every
      t in [T, T+14] (i.e. the correlation has crossed and stays above the
      threshold for 15 consecutive years). If no such T exists, the country
      is reported as "no sustained sync up to 2025".
    - a secondary metric, corr on year-on-year log growth (1-yr differences
      of the logged series), is stored as a robustness check.

Outputs:
    data/wpi_em_sync.json                sync-start table + correlation series
    figures/wpi_em_gold_denominated_sync.jpg   all countries + US, one figure,
                                               sync-start years marked

Usage:
    python scripts/wpi_em_sync.py
"""
import os, sys, json
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, 'data')
FIG = os.path.join(REPO, 'figures')

EM = {'CN': 'China', 'IN': 'India', 'BR': 'Brazil', 'MX': 'Mexico',
      'KE': 'Kenya', 'NG': 'Nigeria', 'EG': 'Egypt', 'ID': 'Indonesia',
      'PK': 'Pakistan', 'MA': 'Morocco'}
CORR_THRESH = 0.5
CORR_WIN = 20          # years in the rolling window
PERSIST = 15           # years the threshold must hold
Y0, Y1 = 1950, 2026    # plot range (US available from 1950 in this dataset)

def load():
    with open(os.path.join(DATA, 'wpi_em_annual.json'), encoding='utf-8') as f:
        em = json.load(f)
    with open(os.path.join(DATA, 'wpi_annual_1800_2026.json'), encoding='utf-8') as f:
        us = json.load(f)
    # US gold-denominated WPI: rebase the 1913=100 series to 2010 = 100.
    g = {int(y): v for y, v in us['annual']['goldwpi']['US'].items() if v is not None}
    base = g.get(2010, g[min(g, key=lambda y: abs(y - 2010))])
    g = {y: 100.0 * v / base for y, v in g.items() if 1950 <= y <= Y1}
    return em, g

def smoothed_log(ys):
    """logged series smoothed with the centred 9-yr MA (pandas convention of
    this repository: min_periods=2 so the window does not shorten the series)."""
    s = pd.Series({y: np.log(v) for y, v in ys.items() if v and v > 0}).sort_index()
    return s.rolling(9, center=True, min_periods=2).mean()

def rolling_corr(x, y):
    """rolling Pearson correlation of two year-indexed series; returns
    {end-year: corr} for every end-year with >=16 matched years in the window."""
    idx = sorted(set(x.index) & set(y.index))
    out = {}
    for t in range(idx[0] + CORR_WIN - 1, idx[-1] + 1):
        w = [yy for yy in idx if t - CORR_WIN + 1 <= yy <= t]
        if len(w) >= 16:
            a, b = x.loc[w].values, y.loc[w].values
            if a.std() > 0 and b.std() > 0:
                out[t] = float(np.corrcoef(a, b)[0, 1])
    return out

def first_sustained(corr, years):
    """earliest t with corr(t) >= threshold for PERSIST consecutive years;
    None if it never holds."""
    above = sorted(t for t in years if corr.get(t, -1) >= CORR_THRESH)
    run = 0
    prev = None
    for t in above:
        run = run + 1 if (prev is not None and t == prev + 1) else 1
        if run >= PERSIST:
            return t - PERSIST + 1
        prev = t
    return None

def main():
    em, us_g = load()
    us_s = smoothed_log(us_g)
    us_growth = us_s.diff()

    out, labels = {}, []
    for cc, nm in EM.items():
        g = {int(y): v for y, v in em['series'][cc]['goldden'].items()}
        s = smoothed_log(g)
        corr_l = rolling_corr(s, us_s)
        corr_g = rolling_corr(s.diff(), us_growth)
        T = first_sustained(corr_l, sorted(corr_l))
        Tg = first_sustained(corr_g, sorted(corr_g))
        out[cc] = {'country': nm, 'span': f'{min(g)}-{max(g)}',
                   'sync_start_level': T, 'sync_start_growth': Tg,
                   'corr_level_end': corr_l.get(max(corr_l)),
                   'corr_growth_end': corr_g.get(max(corr_g)),
                   'corr_level_ts': {str(t): round(v, 3) for t, v in
                                     sorted(corr_l.items())},
                   'corr_growth_ts': {str(t): round(v, 3) for t, v in
                                      sorted(corr_g.items())}}
        labels.append((cc, nm, T))
        print(f'{cc} {nm:10s} span {min(g)}-{max(g)}   '
              f'sync-from (level): {T if T else "none"}   '
              f'(growth): {Tg if Tg else "none"}')

    meta = {'method': f'rolling {CORR_WIN}-y Pearson correlation of centred 9-yr MA of '
                      f'log(gold-denominated index); sync start = first year with corr >= '
                      f'{CORR_THRESH} sustained for {PERSIST} years',
            'us_series': 'US gold-denominated WPI (data/wpi_annual_1800_2026.json, '
                         'rebased 2010=100)',
            'threshold': CORR_THRESH}
    with open(os.path.join(DATA, 'wpi_em_sync.json'), 'w', encoding='utf-8') as f:
        json.dump({'meta': meta, 'countries': out}, f, ensure_ascii=False, indent=1)
    print('wrote data/wpi_em_sync.json')
    draw_figure(em, us_g, labels)

def draw_figure(em, us_g, labels):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors = {'CN': '#e6194B', 'IN': '#f58231', 'BR': '#2ca02c', 'MX': '#ffe119',
              'KE': '#911eb4', 'NG': '#42d4f4', 'EG': '#f032e6', 'ID': '#bfef45',
              'PK': '#3cb44b', 'MA': '#a9a9a9', 'US': '#111111'}
    fig, ax = plt.subplots(figsize=(17, 10.5), dpi=130)
    # thin annual + thick 9-yr MA per country
    ma_handles = []
    for cc, nm in EM.items():
        g = {int(y): v for y, v in em['series'][cc]['goldden'].items()}
        col = colors[cc]
        ax.plot(sorted(g), [g[y] for y in sorted(g)], lw=0.5, color=col, alpha=0.45)
        s = pd.Series(g).sort_index()
        ma = s.rolling(9, center=True, min_periods=2).mean().dropna()
        ax.plot(ma.index, ma.values, lw=2.0, color=col,
                label=f'{cc} - {nm} ({min(g)})')
    # US: black thick MA
    us = pd.Series(us_g).sort_index()
    ma = us.rolling(9, center=True, min_periods=2).mean().dropna()
    ax.plot(sorted(us_g), [us_g[y] for y in sorted(us_g)], lw=0.7, color='#111111', alpha=0.5)
    ax.plot(ma.index, ma.values, lw=3.0, color='#111111', label='US - United States (1950)')
    ax.set_yscale('log')
    ax.set_xlim(Y0 - 1, Y1)
    ax.set_ylim(30, 1500)
    ax.grid(alpha=0.3, which='both')
    ax.set_title('Gold-denominated wholesale prices, PPP-corrected, 10 emerging economies '
                 'vs the United States (2010 = 100, log scale)\n'
                 'solid = 9-yr MA of the annual index; dashed vertical lines mark the first '
                 'year of sustained sync with the US (rolling 20-y corr >= 0.5 for 15 y)',
                 fontsize=12)
    ax.set_xlabel('Year')
    ax.set_ylabel('Index (2010 = 100)')
    # group sync years (several countries share the earliest measurable year)
    groups = {}
    for cc, nm, T in labels:
        if T is None:
            continue
        groups.setdefault(T, []).append(cc)
    for T in sorted(groups):
        ccs = groups[T]
        col = colors[ccs[0]]
        for cc in ccs:
            ax.axvline(T, color=colors[cc], ls='--', lw=1.1, alpha=0.75)
        txt = f'{ccs[0]} {T}' if len(ccs) == 1 else f'{len(ccs)} countries {T}'
        ax.text(T, 1380, txt, fontsize=8.5, fontweight='bold', ha='left',
                bbox=dict(boxstyle='round,pad=0.12', fc='white', alpha=0.85,
                          ec='#888888', lw=0.4), zorder=9)
    ax.text(0.995, 0.012, 'No sync-start marker = never sustainedly correlated '
            '(criterion above)',
            transform=ax.transAxes, ha='right', va='bottom', fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.85, ec='#888888', lw=0.5))
    ax.legend(loc='upper right', fontsize=8.2, ncol=2, framealpha=0.92)
    fig.tight_layout()
    os.makedirs(FIG, exist_ok=True)
    # matplotlib Agg cannot write JPG directly; render PNG via matplotlib and
    # convert to JPG with Pillow.
    import io
    from PIL import Image
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=130)
    plt.close(fig)
    buf.seek(0)
    img = Image.open(buf).convert('RGB')
    img.save(os.path.join(FIG, 'wpi_em_gold_denominated_sync.jpg'), quality=93)
    print('saved figures/wpi_em_gold_denominated_sync.jpg')

if __name__ == '__main__':
    main()
