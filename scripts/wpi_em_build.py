# -*- coding: utf-8 -*-
"""
Build an annual price / PPP / gold-denominated index dataset for ten emerging
economies (China, India, Brazil, Mexico, Kenya, Nigeria, Egypt, Indonesia,
Pakistan, Morocco), 1950s-2025, with the United States kept as the reference
country for the cycle-synchronisation analysis (see wpi_em_sync.py).

Method (one annual series per country):

1. Price index in local currency ("pi"):
   - Brazil: FGV IGP-DI (general price index, ~60% wholesale weights - the
     only long wholesale-based series that could be fetched headlessly) from
     the Banco Central do Brasil SGS API (monthly % changes, Feb 1944+),
     chained into an annual index.
   - All other countries: consumer price index (World Bank, FP.CPI.TOTL,
     2010=100), from the earliest available year (1960 for IN/MX/KE/NG/EG/
     ID/PK/MA; 1986 for China).
   - Caveat: the official WPI/PPI portals of most of these countries (India
     OEA/RBI, Mexico INEGI/Banxico, China NBS, ...) block headless HTTP
     access; a genuine WPI/PPI series is therefore only used for Brazil.
     Swapping in an official WPI/PPI later is a drop-in replacement of the
     "pi" series (the pipeline below is unchanged).

2. PPP conversion factor ("ppp", LCU per international dollar):
   - Official values from the World Bank (PA.NUS.PPP) start in 1990.
   - Before 1990 the factor is backcast with the relative-inflation rule
     ppp(t) = ppp(1990) * [pi(t)/pi(1990)] / [CPI_US(t)/CPI_US(1990)],
     i.e. a constant 1990 real exchange rate (US CPI from FRED CPIAUCSL).
     This is a documented approximation - no official PPP exists before 1990.

3. Gold-denominated, PPP-corrected index ("goldden"):
       goldden(t) = pi(t) / [ppp(t) * gold_usd(t)],  rebased 1913-like: 2010 = 100.
   Because developing currencies often depreciated far beyond purchasing-power
   parity, the PPP factor plays the role of the "fair" exchange rate when
   denomating the index in gold (the raw USD-FX version collapses during
   hyper-depreciation episodes).

Gold: annual average USD/oz 1950-2025 from a monthly market-price CSV
(gold_monthly.csv, default location <repo>/wpi_raw_data/, git-ignored). A
static copy of the derived annual gold series is committed inside
data/wpi_em_annual.json ("gold_usd") so the sync analysis can run without the
raw file.

Usage:
    python scripts/wpi_em_build.py [--gold-csv <path>] [--out data/wpi_em_annual.json]

Output:
    data/wpi_em_annual.json    {"meta": ..., "gold_usd": {...},
                                "us_cpi": {...},
                                "series": {CC: {"pi": {...}, "ppp": {...},
                                                "goldden": {...}}}}
"""
import os, sys, json, argparse, urllib.request, ssl
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(REPO, 'wpi_raw_data')
OUT = os.path.join(REPO, 'data')
os.makedirs(OUT, exist_ok=True)

COUNTRIES = {'CN': 'China', 'IN': 'India', 'BR': 'Brazil', 'MX': 'Mexico',
             'KE': 'Kenya', 'NG': 'Nigeria', 'EG': 'Egypt', 'ID': 'Indonesia',
             'PK': 'Pakistan', 'MA': 'Morocco'}
BASE_YEAR = 2010

def fetch(url, timeout=60, verify=True):
    ctx = ssl.create_default_context() if verify else ssl._create_unverified_context()
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
        return r.read().decode('utf-8', 'ignore')

def wb_indicator(cc, ind, y0, y1):
    """World Bank indicator {year: value}; empty dict when the indicator or
    range is unavailable."""
    url = (f'https://api.worldbank.org/v2/country/{cc}/indicator/{ind}'
           f'?format=json&date={y0}:{y1}&per_page=400')
    t = fetch(url)
    d = json.loads(t)
    rows = d[1] if isinstance(d, list) and len(d) > 1 else []
    out = {}
    for r in rows:
        if r.get('value') is not None:
            out[int(r['date'])] = float(r['value'])
    return out

def wb_cpi(cc):
    return wb_indicator(cc, 'FP.CPI.TOTL', 1950, 2025)

def wb_ppp(cc):
    return wb_indicator(cc, 'PA.NUS.PPP', 1990, 2025)

def fred_annual_us_cpi():
    """US CPI-U annual averages (FRED CPIAUCSL, monthly 1947-2026)."""
    t = fetch('https://fred.stlouisfed.org/graph/fredgraph.csv?id=CPIAUCSL', timeout=40)
    out = {}
    for ln in t.strip().splitlines()[1:]:
        d, v = ln.split(',')
        if not v.strip() or v == '.':
            continue
        y = int(d[:4])
        out.setdefault(y, []).append(float(v))
    return {y: float(np.mean(v)) for y, v in out.items()}

def brazil_igp_annual():
    """FGV IGP-DI: chain the monthly % changes (BCB SGS series 190, Feb 1944+)
    into an annual index (Feb 1944 = 100)."""
    t = fetch('https://api.bcb.gov.br/dados/serie/bcdata.sgs.190/dados'
              '?formato=json&dataInicial=01/01/1944&dataFinal=31/12/2025')
    rows = json.loads(t)
    monthly = {}          # (year, month) -> change in %
    for r in rows:
        d, v = r['data'], float(r['valor'])
        mm, yy = int(d[3:5]), int(d[6:10])
        monthly[(yy, mm)] = v
    idx = {}
    level = 100.0
    for (yy, mm) in sorted(monthly):
        level *= 1.0 + monthly[(yy, mm)] / 100.0
        idx.setdefault(yy, []).append(level)
    return {y: float(np.mean(v)) for y, v in idx.items()}

def gold_annual(gold_csv):
    if not gold_csv or not os.path.exists(gold_csv):
        raise SystemExit(f'gold CSV not found: {gold_csv} - provide --gold-csv '
                         '(expected: date,gold monthly, from the commodity project)')
    t = open(gold_csv, encoding='utf-8').read()
    out = {}
    for ln in t.strip().splitlines()[1:]:
        d, v = ln.split(',')
        if not v.strip():
            continue
        out.setdefault(int(d[:4]), []).append(float(v))
    return {y: float(np.mean(v)) for y, v in out.items() if y <= 2025}

def main():
    ap = argparse.ArgumentParser(description='Build EM annual price/PPP/gold series.')
    ap.add_argument('--gold-csv', default=os.path.join(RAW, 'gold_monthly.csv'))
    ap.add_argument('--out', default=os.path.join(OUT, 'wpi_em_annual.json'))
    args = ap.parse_args()

    us_cpi = fred_annual_us_cpi()
    gold = gold_annual(args.gold_csv)
    yrs_g = sorted(gold)

    print('fetching World Bank CPI / PPP and Brazil IGP-DI ...')
    cpi, ppp = {}, {}
    for cc in COUNTRIES:
        cpi[cc] = wb_cpi(cc)
        ppp[cc] = wb_ppp(cc)
        print(f'  {COUNTRIES[cc]:9s} CPI {min(cpi[cc]):.0f}-{max(cpi[cc]):.0f} '
              f'({len(cpi[cc])}y)   PPP {min(ppp[cc]):.0f}-{max(ppp[cc]):.0f} ({len(ppp[cc])}y)')
    br = brazil_igp_annual()
    print(f'  Brazil   IGP-DI (FGV wholesale-based) {min(br):.0f}-{max(br):.0f} ({len(br)}y)')

    # ---- PPP backcast (pre-1990) with the relative-inflation rule ----
    for cc in COUNTRIES:
        pi = br if cc == 'BR' else cpi[cc]
        anchor = 1990
        for y in range(min(pi), 1990):
            if y in pi and anchor in pi and y in us_cpi and anchor in us_cpi:
                ppp[cc][y] = (ppp[cc][anchor] * (pi[y] / pi[anchor])
                              / (us_cpi[y] / us_cpi[anchor]))
        ppp[cc] = {y: ppp[cc][y] for y in sorted(ppp[cc])}

    # ---- gold-denominated, PPP-corrected index (2010 = 100) ----
    series = {}
    for cc in COUNTRIES:
        pi = br if cc == 'BR' else cpi[cc]
        gd = {}
        for y in sorted(set(pi) & set(ppp[cc]) & set(gold)):
            gd[y] = pi[y] / (ppp[cc][y] * gold[y])
        base = gd.get(BASE_YEAR)
        if not base:
            base = gd.get(min(gd, key=lambda y: abs(y - BASE_YEAR)))
        gd = {y: 100.0 * v / base for y, v in gd.items()}
        series[cc] = {'pi': pi, 'ppp': ppp[cc], 'goldden': gd}
        print(f'  {COUNTRIES[cc]:9s} goldden {min(gd):.0f}-{max(gd):.0f} '
              f'({len(gd)}y, base {BASE_YEAR})')

    meta = {
        'title': 'Annual local-currency price index, PPP conversion factor and '
                 'PPP-corrected gold-denominated index, 10 emerging economies + US reference',
        'base_year': BASE_YEAR,
        'goldden_formula': 'goldden(t) = pi(t) / [ppp(t) * gold_usd(t)], rebased to 2010 = 100',
        'price_index': {
            'BR': 'FGV IGP-DI chained from monthly % changes (BCB SGS s.190, 1944+): '
                  'general price index with ~60% wholesale weights (WPI proxy)',
            'default': 'consumer price index, World Bank FP.CPI.TOTL (2010=100). Official '
                       'WPI/PPI portals (India OEA/RBI, Mexico INEGI/Banxico, China NBS, '
                       'etc.) are not reachable headlessly; swap-in possible via the pi series'},
        'coverage': {cc: f'{min(series[cc]["goldden"])}-{max(series[cc]["goldden"])}'
                     for cc in COUNTRIES},
        'ppp_method': 'World Bank PA.NUS.PPP (LCU per international dollar) from 1990; '
                      'pre-1990 backcast: ppp(t) = ppp(1990) * [pi(t)/pi(1990)] / '
                      '[CPI_US(t)/CPI_US(1990)] (constant 1990 real exchange rate; US CPI '
                      'from FRED CPIAUCSL). No official PPP exists before 1990.',
        'gold': 'annual average USD/oz, 1950-2025, from monthly market prices '
                '(gold_monthly.csv; London fixing / COMEX; cached from the commodity project)',
        'us_reference': 'United States: WPI in gold terms from data/wpi_annual_1800_2026.json '
                        '(rebased to 2010=100 in wpi_em_sync.py)',
    }
    with open(args.out, 'w', encoding='utf-8') as f:
        json.dump({'meta': meta, 'gold_usd': {str(y): v for y, v in gold.items()},
                   'us_cpi': {str(y): v for y, v in us_cpi.items()},
                   'series': {cc: {k: {str(y): v for y, v in d.items()}
                                   for k, d in series[cc].items()}
                              for cc in series}},
                  f, ensure_ascii=False, indent=1)
    print('wrote', args.out)

if __name__ == '__main__':
    main()
