# -*- coding: utf-8 -*-
"""
Plot the local-currency WPI trends, 1800-2026, one subplot per country, with
centred 3/6/9-yr moving averages.

- Every country series is rebased so that its own 1913 level = 100 (log y).
- The pre-1940 span is highlighted by a shaded background; peaks and troughs
  of the centred 9-yr MA (min. 3-yr spacing, swing >= 12%, alternating) are
  labelled there only.
- Germany: the 1923 hyperinflation (annual average 1.26e11) would flatten all
  other German history on a full-range log axis, so the German panel caps the
  y-axis at 1e5 and the off-scale peak is annotated with an arrow and a note.

Input:
    data/wpi_annual_1800_2026.json    annual local-currency WPI

Outputs:
    figures/wpi_chart_local_currency.png    one figure, 3x2 subplot grid
    data/wpi_pre1940_cycles_local.json      pre-1940 peak/trough table

Usage:
    python scripts/wpi_plot_local_currency.py
"""
import os, sys, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
DATA = os.path.join(REPO, 'data')
FIG = os.path.join(REPO, 'figures')
os.makedirs(FIG, exist_ok=True)

plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei']
plt.rcParams['axes.unicode_minus'] = False

NAMES = {'US': 'United States', 'UK': 'United Kingdom', 'DE': 'Germany',
         'JP': 'Japan', 'FR': 'France'}
ORDER = ['US', 'UK', 'DE', 'JP', 'FR']
PEAK_C, TROUGH_C = '#d62728', '#1f77b4'
MA_STYLE = [(3, '#1f77b4', '3-yr MA', 1.5), (6, '#ff7f0e', '6-yr MA', 1.5),
            (9, '#d62728', '9-yr MA', 2.4)]
DE_CAP = 1e5   # Germany panel y-axis cap (1923 hyperinflation peak is off scale)

def moving_avg(s, w):
    return s.rolling(w, center=True, min_periods=2).mean()

def find_cycles(s):
    """Same detection as scripts/wpi_build_series.py: local extrema of the
    centred 9-yr MA, min. 3-yr spacing, swing >= 12%, alternating; a first
    trough within 5 years of the series start is discarded as an edge artefact.
    Returns [(year, type, ma9-value)]."""
    m9 = s.rolling(9, center=True, min_periods=9).mean().dropna()
    raw = []
    for y in m9.index:
        window = m9[(m9.index >= y - 2) & (m9.index <= y + 2)]
        if len(window) < 5:
            continue
        v = float(m9[y])
        is_peak = v >= window.max()
        is_trough = v <= window.min()
        if not (is_peak or is_trough):
            continue
        typ = 'peak' if is_peak else 'trough'
        if raw and (y - raw[-1][0]) < 3:
            if typ == raw[-1][1] and ((typ == 'peak' and v > raw[-1][2]) or
                                      (typ == 'trough' and v < raw[-1][2])):
                raw[-1] = (y, typ, v)
            continue
        raw.append((y, typ, v))
    start = int(m9.index[0])
    pts = []
    for (y, typ, v) in raw:
        if not pts:
            if typ == 'trough' and y - start < 5:
                continue
            pts.append((y, typ, v)); continue
        last_y, last_t, last_v = pts[-1]
        if typ == last_t:
            if (typ == 'peak' and v > last_v) or (typ == 'trough' and v < last_v):
                pts[-1] = (y, typ, v)
            continue
        if typ == 'peak' and v <= last_v: continue
        if typ == 'trough' and v >= last_v: continue
        if abs(v - last_v) / max(abs(last_v), 1e-9) < 0.12: continue
        pts.append((y, typ, round(v, 1)))
    if not pts and len(m9) > 3:
        pts.append((int(m9.index[0]), 'start', round(float(m9.iloc[0]), 1)))
    return pts

def main():
    with open(os.path.join(DATA, 'wpi_annual_1800_2026.json'), encoding='utf-8') as f:
        ann_j = json.load(f)
    fig, axs = plt.subplots(3, 2, figsize=(16, 12), dpi=110)
    fig.subplots_adjust(hspace=0.22, wspace=0.16, top=0.955, bottom=0.05,
                        left=0.055, right=0.985)
    print('%-4s  pre-1940 cycles (year P/T on the 9-yr MA)' % 'CCY')
    cyc_save = {}

    for idx, c in enumerate(ORDER):
        ax = axs[idx // 2][idx % 2]
        # Rebase the published local-currency index to its own 1913 = 100.
        yr = [int(y) for y in ann_j['annual']['wpi'][c]]
        val = ann_j['annual']['wpi'][c]
        d = pd.DataFrame({'year': yr,
                          'wpi': [val[str(y)] if val.get(str(y)) is not None else np.nan
                                  for y in yr]}).dropna().reset_index(drop=True)
        base13 = float(d[d.year == 1913]['wpi'].iloc[0])
        s = d.set_index('year')['wpi'] / base13 * 100.0
        y0, y1 = int(s.index.min()), int(s.index.max())

        cyc = [(y, t, v) for (y, t, v) in find_cycles(s) if y <= 1939]
        cyc_save[c] = cyc
        m9 = s.rolling(9, center=True, min_periods=9).mean().dropna()

        # ---- draw ----
        ax.axvspan(max(1800, y0), 1940, color='#fff3bf', alpha=0.35, lw=0, zorder=0)
        ax.axvline(1940, color='#999999', ls=':', lw=1.0, zorder=1)
        ax.plot(s.index, s.values, color='#bbbbbb', lw=0.9, zorder=2)
        for w, col, lab, lw in MA_STYLE:
            ma = moving_avg(s, w)
            ax.plot(ma.index, ma.values, color=col, lw=lw, zorder=3)

        ax.set_yscale('log')
        if c == 'DE':
            ax.set_ylim(10 ** 1.7, DE_CAP)
        else:
            ax.set_xlim(y0 - 1, 2026.5)
            ax.relim(); ax.autoscale_view()
        ax.autoscale(False)
        ax.set_xlim(y0 - 1, 2026.5)

        ylo, yhi = ax.get_ylim()
        span = np.log10(yhi) - np.log10(ylo)
        for (yy, typ, v0) in cyc:
            if yy not in m9.index:
                continue
            v = float(v0)
            pos = (np.log10(v) - np.log10(ylo)) / span
            if c == 'DE' and v > DE_CAP:
                # Off-scale peak (9-yr MA ~1.4e10 over 1919-27): draw the marker
                # at the axis top edge and explain in the panel note below.
                ax.plot(yy, 10 ** (np.log10(DE_CAP) - 0.06), marker='^',
                        color=PEAK_C, ms=11, zorder=6)
                ax.annotate(f'{yy} (off scale)',
                            xy=(yy, 10 ** (np.log10(DE_CAP) - 0.06)),
                            xytext=(6, 3), textcoords='offset points', ha='left',
                            va='bottom', fontsize=8.5, fontweight='bold',
                            color=PEAK_C,
                            bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.9,
                                      ec='none'), zorder=8)
                continue
            col = PEAK_C if typ == 'peak' else TROUGH_C
            m = '^' if typ == 'peak' else 'v'
            ax.plot(yy, v, marker=m, color=col, ms=11, zorder=6)
            if typ == 'peak':
                dy = 9 if pos < 0.84 else -14     # peak labels usually above
            else:
                dy = -14 if pos > 0.16 else 9     # trough labels usually below
            ax.annotate(f'{yy}', xy=(yy, v), xytext=(0, dy), textcoords='offset points',
                        ha='center', fontsize=9.5, fontweight='bold', color=col,
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.9,
                                  ec='none'), clip_on=False, zorder=8)

        ax.set_title(f'{NAMES[c]}  |  local-currency WPI, {y0}-{y1}  (1913 = 100)',
                     fontsize=12)
        ax.grid(alpha=0.3, which='both')
        if c == 'DE':
            ax.text(0.012, 0.985, '1923 annual avg = 1.26e11 (Dec peak 4.2e12)\n'
                    'above the 1e5 axis cap;  9-yr MA ~1.4e10 over 1919-27.\n'
                    '1924 Rentenmark reform ~137 (pre-war level).',
                    transform=ax.transAxes, va='top', ha='left', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.35', fc='white', alpha=0.92,
                              ec='#888888', lw=0.5), zorder=9)
            ax.text(0.012, 0.02, '1914-18 WWI inflation  |  1945-47 no data  |'
                    '  1948+ West Germany', transform=ax.transAxes, va='bottom',
                    ha='left', fontsize=7.6, color='#444444', zorder=9)
        if idx % 2 == 0:
            ax.set_ylabel('Index (1913 = 100)', fontsize=9)
        if idx // 2 == 2:
            ax.set_xlabel('Year', fontsize=10)
        else:
            ax.tick_params(labelbottom=False)
        print('%-4s  %s' % (c, '  '.join(f'{y}{"P" if t=="peak" else "T"}({v:.0f})'
                                         for (y, t, v) in cyc)))

    # ---- 6th panel: legend and notes ----
    ax = axs[2][1]
    ax.axis('off')
    proxies = [plt.Line2D([], [], color='#bbbbbb', lw=1.0,
                          label='Annual WPI (local currency)')]
    for w, col, lab, lw in MA_STYLE:
        proxies.append(plt.Line2D([], [], color=col, lw=lw, label=lab))
    proxies += [plt.Line2D([], [], marker='^', ls='', color=PEAK_C, ms=10,
                           label='Peak (9-yr MA max)'),
                plt.Line2D([], [], marker='v', ls='', color=TROUGH_C, ms=10,
                           label='Trough (9-yr MA min)')]
    ax.legend(handles=proxies, loc='upper center', fontsize=10, framealpha=0.9,
              ncol=2, columnspacing=1.2, handlelength=1.6)
    notes = [
        'Each panel: wholesale price index in the country\'s own currency,\n'
        'rebased to 1913 = 100.  Log scale.',
        'Coverage: US/UK 1800-2026, Germany 1851-2026 (no data 1945-47),\n'
        'Japan 1887-2026, France 1913-2025.',
        'Shaded zone 1800-1940: peaks and troughs of the centred 9-yr MA\n'
        'are labelled there only (min. 3-yr spacing, swing >= 12%, alternating).',
        'Germany: hyperinflation peak 1923 (annual avg 1.26e11, i.e. 126 billion\n'
        'times the 1913 level) lies above the axis cap of 1e5; see DE panel note.',
        'Note: historical WPI series are spliced across sources and definitions\n'
        '(wholesale -> producer prices); level breaks of a few percent possible.',
    ]
    for i, t in enumerate(notes):
        ax.text(0.03, 0.66 - i * 0.135, t, transform=ax.transAxes, fontsize=9.5,
                va='top', ha='left')

    fig.suptitle('Wholesale Price Index in Local Currency, 1800-2026'
                 '  (own currency, 1913 = 100)', fontsize=15)

    # ---- label-box overlap check ----
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    n_warn = 0
    for ax in axs.flat:
        if not ax.axison:
            continue
        labs = [t for t in ax.texts if t.get_text().strip()]
        boxes = [(t, t.get_window_extent(rend)) for t in labs]
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i][1], boxes[j][1]
                if not (a.x1 < b.x0 or b.x1 < a.x0 or a.y1 < b.y0 or b.y1 < a.y0):
                    print(f'  WARN label overlap: "{boxes[i][0].get_text()}"'
                          f' vs "{boxes[j][0].get_text()}"')
                    n_warn += 1
    print('label overlap check: %d warning(s)' % n_warn)

    fig.savefig(os.path.join(FIG, 'wpi_chart_local_currency.png'), bbox_inches='tight')
    plt.close(fig)
    print('saved figures/wpi_chart_local_currency.png')

    rows = []
    for c in ORDER:
        for (y, t, v) in cyc_save[c]:
            rows.append({'country': c, 'year': y, 'type': t, 'ma9_value': v})
    with open(os.path.join(DATA, 'wpi_pre1940_cycles_local.json'), 'w', encoding='utf-8') as f:
        json.dump({'title': 'Pre-1940 peak/trough cycles of the local-currency WPI '
                            '(labelled on figures/wpi_chart_local_currency.png)',
                   'detection': 'same 9-yr MA method as wpi_cycles_gold.json; only '
                                'extrema with year < 1940 are listed',
                   'cycles': rows}, f, ensure_ascii=False, indent=1)
    print('saved data/wpi_pre1940_cycles_local.json')

if __name__ == '__main__':
    main()
