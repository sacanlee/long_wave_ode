# -*- coding: utf-8 -*-
"""
Build the gold-denominated wholesale price index (WPI) series:
US / UK / Germany / Japan / France, 1800-2026.

- 1800-1981: annual; 1982+: quarterly (monthly PPI / gold price / FX aggregated).
- Gold-denominated index = WPI(local currency) / [FX(local currency per USD)
  x gold price(USD/oz)] = commodity price / gold price in local currency.
  Base year: 1913 = 100.
- Peak/trough cycle detection on the centred 9-yr moving average of the annual
  series (see data/wpi_cycles_gold.json).

Usage:
    python scripts/wpi_build_series.py [--raw-dir <dir>]

--raw-dir points to a directory with the externally downloaded source files
(raw FRED CSVs, BOE/INSEE/Destatis workbooks, BOJ PDFs, etc.). The file list
and download locations are documented in README.md (section "WPI input data").
If omitted, the default is <repo>/wpi_raw_data/.

Outputs (regenerated results, repo JSON convention):
    data/wpi_annual_1800_2026.json     annual wpi / fx / goldwpi / gold price
    data/wpi_quarterly_1982_2026.json  quarterly gold-denominated index
    data/wpi_cycles_gold.json          peak/trough cycles on the 9-yr MA
"""
import os, re, sys, json, argparse
import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # repo root
RAW = os.path.join(REPO, 'wpi_raw_data')     # external downloads (not in git)
OUT = os.path.join(REPO, 'data')             # regenerated JSON results
os.makedirs(OUT, exist_ok=True)

YEARS = list(range(1800, 2027))

def read_fred_csv(fn):
    df = pd.read_csv(os.path.join(RAW, fn))
    df.columns = ['date', 'value']
    df['date'] = pd.to_datetime(df['date'])
    return df.set_index('date')['value'].dropna().sort_index()

# ============ 1. Gold price USD/oz ============
# Civil-war greenback premium over gold (Mitchell 1908 / HSUS), % of par.
GOLD_PREMIUM = {1861: 100.0, 1862: 113.2, 1863: 145.3, 1864: 203.3, 1865: 157.3,
                1866: 141.6, 1867: 138.2, 1868: 139.8, 1869: 133.6, 1870: 114.5,
                1871: 112.9, 1872: 112.3, 1873: 113.9, 1874: 112.4, 1875: 116.8,
                1876: 111.5, 1877: 105.3, 1878: 101.7, 1879: 100.0}

def gold_annual_series():
    """Annual USD/oz gold price: official prices before 1950, monthly market
    averages from the cached gold_monthly.csv (1950-2026)."""
    g = {}
    for y in YEARS:
        if y <= 1833: g[y] = 19.39
        elif y <= 1861: g[y] = 20.67
        elif y <= 1878: g[y] = 20.67 * GOLD_PREMIUM[y] / 100.0
        elif y <= 1933: g[y] = 20.67
        elif y <= 1967: g[y] = 35.0
        else: g[y] = np.nan
    gold_m = pd.read_csv(os.path.join(RAW, 'gold_monthly.csv'), parse_dates=['date'])
    gold_m = gold_m.set_index('date')['gold'].dropna()
    for y in range(1950, 2027):
        sel = gold_m[gold_m.index.year == y]
        if len(sel): g[y] = float(sel.mean())
    return g, gold_m

# ============ 2. Exchange rates (local currency per USD) ============
def parse_mw(fn):
    """Parse a MeasuringWorth currency-comparison HTML table into {year: rate}."""
    html = open(os.path.join(RAW, fn), encoding='utf-8', errors='ignore').read()
    rows = re.findall(r'<tr[^>]*>(.*?)</tr>', html, re.S)
    out = {}
    for r in rows:
        cells = re.findall(r'<t[dh][^>]*>(.*?)</t[dh]>', r, re.S)
        cells = [re.sub(r'<[^>]+>', '', x).strip() for x in cells]
        if len(cells) == 2 and cells[0].isdigit():
            v = cells[1].split(' ')[0].replace(',', '')
            try: out[int(cells[0])] = float(v)
            except Exception: pass
    return out

def parse_ltes_fx():
    """Hitotsubashi LTES workbook, table 26: USD per 100 yen, 1874-1941,
    converted to yen per USD."""
    wb = __import__('openpyxl').load_workbook(os.path.join(RAW, 'ltes14.xlsx'),
                                              read_only=True, data_only=True)
    ws = wb['第26表']
    out = {}
    for r in ws.iter_rows(values_only=True):
        if r[3] is None: continue
        y = str(r[3]).strip()
        if not y.isdigit(): continue
        v = r[4]
        if isinstance(v, (int, float)) and v > 0:
            out[int(y)] = 100.0 / float(v)
    wb.close()
    return out

def fx_annual_series():
    fx = {}
    mw_de = parse_mw('mw_Germany.html')
    mw_fr = parse_mw('mw_France.html')
    mw_jp = parse_mw('mw_Japan.html')
    ltes = parse_ltes_fx()

    # UK: BOE A33 USD/GBP (1791-2016) -> GBP/USD; 2017+ FRED EXUSUK.
    uk = {}
    wb = __import__('openpyxl').load_workbook(os.path.join(RAW, 'boe_millennium.xlsx'),
                                              read_only=True, data_only=True)
    ws = wb['A33. Exchange rate data']
    for r in ws.iter_rows(values_only=True):
        y, v = r[0], r[1]
        if isinstance(y, (int, float)) and isinstance(v, (int, float)):
            uk[int(y)] = 1.0 / float(v)
    wb.close()
    exusuk = read_fred_csv('fred_EXUSUK.csv')
    for y in range(2017, 2027):
        sel = exusuk[exusuk.index.year == y]
        if len(sel): uk[y] = float(sel.mean())

    # DE: <=1912 gold parity 4.198; 1923 hyperinflation annual average
    # (Dec 1923 peak 4.2e12 marks per USD); 1924 Rentenmark parity 4.2;
    # 1913-1970 MeasuringWorth; 1971+ FRED EXGEUS/DEXUSEU.
    de = {}
    for y in YEARS:
        if y <= 1912: de[y] = 4.198
        elif y == 1923: de[y] = 4.2e11
        elif y == 1924: de[y] = 4.2
        elif y <= 1940: de[y] = mw_de.get(y, np.nan)
        elif y <= 1945: de[y] = 2.5
        elif y <= 1947: de[y] = 10.0
        elif y == 1948: de[y] = 3.33
        elif y == 1949: de[y] = 3.62
        elif y <= 1970: de[y] = mw_de.get(y, np.nan)
    exgeus = read_fred_csv('fred_EXGEUS.csv')
    dexeu = read_fred_csv('fred_DEXUSEU.csv')
    de_m = pd.concat([exgeus, dexeu])
    de_m = de_m[~de_m.index.duplicated(keep='first')].sort_index()
    for y in range(1971, 2027):
        sel = de_m[de_m.index.year == y]
        if len(sel): de[y] = float(sel.mean())

    # JP: LTES 1874-1941; official rates 1942-48 (4.25/15/50/270); 1949-1970
    # MeasuringWorth (~360); 1971+ FRED EXJPUS. LTES lacks 1941 (use official
    # 4.27) and 1923 (Kanto earthquake; gap filled by linear interpolation).
    jp = {}
    for y in YEARS:
        if y <= 1940: jp[y] = ltes.get(y, np.nan)
        elif y <= 1945: jp[y] = 4.25
        elif y == 1946: jp[y] = 15.0
        elif y == 1947: jp[y] = 50.0
        elif y == 1948: jp[y] = 270.0
        elif y <= 1970: jp[y] = mw_jp.get(y, 360.0)
    if 1941 not in jp: jp[1941] = 4.27
    yrs = sorted(y for y in jp if y <= 1941 and not np.isnan(jp[y]))
    for y in range(1887, 1942):
        if y not in jp or np.isnan(jp[y]):
            lo = max([yy for yy in yrs if yy < y], default=None)
            hi = min([yy for yy in yrs if yy > y], default=None)
            if lo and hi:
                jp[y] = jp[lo] + (jp[hi] - jp[lo]) * (y - lo) / (hi - lo)
    exjpus = read_fred_csv('fred_EXJPUS.csv')
    for y in range(1971, 2027):
        sel = exjpus[exjpus.index.year == y]
        if len(sel): jp[y] = float(sel.mean())

    # FR in old-franc terms (same scale as the WPI): 1960 redenomination x100;
    # 2002 euro x655.957 old francs per euro.
    fr = {}
    for y in YEARS:
        if y <= 1912: fr[y] = 5.1826
        elif y <= 1940: fr[y] = mw_fr.get(y, np.nan)
        elif y <= 1944: fr[y] = 50.0
        elif y <= 1959: fr[y] = mw_fr.get(y, np.nan)
        elif y <= 2001: fr[y] = mw_fr.get(y, np.nan) * 100.0 if y in mw_fr else np.nan
        elif y <= 2026: fr[y] = dexeu.get(pd.Timestamp(y, 12, 1), np.nan) * 655.957 if len(dexeu[dexeu.index.year == y]) else np.nan
    exfrus = read_fred_csv('fred_EXFRUS.csv') * 100.0   # new francs -> old francs
    fr_m = pd.concat([exfrus, dexeu * 655.957])
    fr_m = fr_m[~fr_m.index.duplicated(keep='first')].sort_index()
    for y in range(1971, 2027):
        sel = fr_m[fr_m.index.year == y]
        if len(sel): fr[y] = float(sel.mean())

    fx = {'US': {y: 1.0 for y in YEARS}, 'UK': uk, 'DE': de, 'JP': jp, 'FR': fr}
    fx_m = {'UK': exusuk, 'DE': de_m, 'JP': exjpus, 'FR': fr_m,
            'US': pd.Series(1.0, index=pd.date_range('1971-01-01', '2026-12-01', freq='MS'))}
    return fx, fx_m

# ============ 3. WPI (local currency) ============
# US wholesale prices 1800-1850, Warren & Pearson (HSUS E23-42); values are
# near a 1913=100 scale and are re-aligned to the 1850+ segment below.
US_WP = {1800:129,1801:141,1802:118,1803:120,1804:126,1805:139,1806:132,1807:129,1808:134,1809:131,
1810:131,1811:128,1812:135,1813:158,1814:182,1815:165,1816:139,1817:135,1818:143,1819:130,
1820:106,1821:98,1822:102,1823:99,1824:95,1825:100,1826:94,1827:93,1828:91,1829:92,
1830:91,1831:92,1832:89,1833:91,1834:90,1835:98,1836:105,1837:104,1838:104,1839:112,
1840:95,1841:89,1842:81,1843:76,1844:79,1845:84,1846:82,1847:92,1848:82,1849:82,1850:84}

# German wholesale prices 1914-1944, Statistisches Reichsamt (1913=100; via
# FRB bulletins and the literature), including the 1922-23 hyperinflation.
REICHSAMT = {1914:125,1915:148,1916:151,1917:203,1918:245,1919:803,1920:1440,1921:3487,
             1922:147500,1923:1.26e11,1924:137,1925:142,1926:134,1927:138,1928:140,
             1929:137,1930:125,1931:111,1932:97,1933:93,1934:98,1935:102,1936:104,
             1937:106,1938:106,1939:107,1940:111,1941:112,1942:114,1943:116,1944:116}
# West Germany 1948-1957 (1913=100; literature approximation, Bundesbank
# long series; +/-10% uncertainty).
DE_1948_57 = {1948: 117, 1949: 118, 1950: 122, 1951: 145, 1952: 145, 1953: 139, 1954: 136,
              1955: 140, 1956: 141, 1957: 142}
# Japan 1945-1959, BOJ/MOF national WPI (1934-36 = 100 base).
JP_BOJ = {1945: 400, 1946: 1600, 1947: 4800, 1948: 12800, 1949: 20900, 1950: 24700, 1951: 34300,
          1952: 34900, 1953: 35200, 1954: 34300, 1955: 34300, 1956: 35800, 1957: 36800, 1958: 34400, 1959: 34800}
# France 1940-45 (FRB bulletin, 1938=100) and 1946-52 (INSEE, 1938=100).
FRB_FR = {1940: 139, 1941: 171, 1942: 201, 1943: 234, 1944: 265, 1945: 375}
INSEE_FR = {1946: 606, 1947: 819, 1948: 1711, 1949: 2114, 1950: 2289, 1951: 2924, 1952: 3063}

def parse_imes(fn):
    """BOJ IMES historical Tokyo wholesale-price monthly CSVs (Shift-JIS)."""
    out = {}
    with open(os.path.join(RAW, fn), encoding='shift_jis', errors='replace') as f:
        for ln in f:
            p = ln.rstrip('\n').split(',')
            d = p[0].strip()
            if len(d) >= 6 and d[:4].isdigit() and '.' in d[:6]:
                try:
                    y, m = d.split('.')
                    y, m = int(y), int(m)
                    if len(p) > 1 and p[1].strip():
                        out[(y, m)] = float(p[1])
                except Exception:
                    pass
    return out

def read_ons_gb7s():
    """ONS PPI dataset CSV: row 1 = title, row 2 = CDID headers, data from row
    4; returns the GB7S (output producer prices) monthly series."""
    import csv as _csv
    with open(os.path.join(RAW, 'ons_ppi.csv'), encoding='utf-8') as f:
        rows = list(_csv.reader(f))
    cdid = rows[1]
    idx = cdid.index('GB7S')
    MON = {'JAN':1,'FEB':2,'MAR':3,'APR':4,'MAY':5,'JUN':6,'JUL':7,'AUG':8,'SEP':9,'OCT':10,'NOV':11,'DEC':12}
    out = {}
    for r in rows[3:]:
        d, v = r[0], r[idx]
        if not d or not v: continue
        m = re.match(r'(\d{4}) (JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)', d)
        if m:
            out[pd.Timestamp(int(m.group(1)), MON[m.group(2)], 1)] = float(v)
    s = pd.Series(out)
    s.index = pd.DatetimeIndex(s.index)
    return s.sort_index()

def read_destatis():
    """Destatis table 61241-b01: monthly 'Gewerbliche Erzeugnisse insgesamt'
    (2021=100), 2010-2026. The Destatis server is very slow; fetch the xlsx
    with scripts/wpi_dl_destatis.py (plain curl downloads get truncated)."""
    wb = __import__('openpyxl').load_workbook(os.path.join(RAW, 'destatis_ppi3.xlsx'),
                                              read_only=True, data_only=True)
    ws = wb['61241-b01']
    rows = list(ws.iter_rows(values_only=True))
    out = {}
    MON = {'Januar':1,'Februar':2,'März':3,'April':4,'Mai':5,'Juni':6,'Juli':7,'August':8,'September':9,'Oktober':10,'November':11,'Dezember':12}
    for r in rows[4:]:
        d = r[0]
        if not isinstance(d, str) or d == 'Ende der Tabelle': continue
        m = re.match(r'(\w+) (\d{4})', d)
        if not m or m.group(1) not in MON: continue
        y = int(m.group(2))
        if r[1] is not None and isinstance(r[1], (int, float)):
            out[(y, MON[m.group(1)])] = float(r[1])
    wb.close()
    return out

def read_insee_ippi():
    """INSEE IPPI monthly xlsx (2015=100), 2005-2025-10."""
    from datetime import datetime
    wb = __import__('openpyxl').load_workbook(os.path.join(RAW, 'insee_ippi.xlsx'),
                                              read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    out = {}
    for r in ws.iter_rows(values_only=True):
        d = r[0]
        if isinstance(d, datetime) and isinstance(r[1], (int, float)):
            out[(d.year, d.month)] = float(r[1])
    wb.close()
    return out

def parse_boj_tail():
    """BOJ CGPI month-on-month change tables in the monthly-report PDFs
    (2024-2026); returns {(year, month): change in %} for the newest months."""
    out = {}
    MONTHS = {'Jan.':1,'Feb.':2,'Mar.':3,'Apr.':4,'May':5,'June':6,'July':7,'Aug.':8,'Sep.':9,'Oct.':10,'Nov.':11,'Dec.':12}
    for fn in ['cgpi2404.pdf', 'cgpi2504.pdf', 'cgpi2607.pdf']:
        path = os.path.join(RAW, fn)
        if not os.path.exists(path): continue
        doc = __import__('pymupdf').open(path)
        t = doc[0].get_text()
        tokens = t.split()
        cur_year = None
        i = 0
        while i < len(tokens):
            tok = tokens[i]
            if re.match(r'^\d{4}/$', tok):
                cur_year = int(tok[:4]); i += 1; continue
            if tok in MONTHS and cur_year:
                mm = MONTHS[tok]
                nums = []
                j = i + 1
                while j < len(tokens):
                    t2 = tokens[j]
                    if t2 in MONTHS or re.match(r'^\d{4}/$', t2): break
                    try:
                        nums.append(float(t2))
                    except ValueError:
                        pass
                    j += 1
                if nums:
                    out[(cur_year, mm)] = nums[0]
                i = j; continue
            i += 1
        doc.close()
    return out

def wpi_annual_and_monthly():
    """Return (annual {country: {year: wpi}}, monthly {country: pd.Series}).
    Local-currency index on each series' published base (UK stays 2015=100)."""
    annual, monthly = {}, {}

    # ===== US =====
    a = read_fred_csv('fred_M0448AUSM323NNBR.csv')   # 1850-1894
    b = read_fred_csv('fred_M0448BUSM336NNBR.csv')   # 1890-1914
    p = read_fred_csv('fred_PPIACO.csv')             # 1913-2026 (1982=100)
    r_ab = (b[(b.index >= '1890-01-01') & (b.index <= '1894-12-01')] /
            a[(a.index >= '1890-01-01') & (a.index <= '1894-12-01')]).mean()
    r_bp = (p[(p.index >= '1913-01-01') & (p.index <= '1914-12-01')] /
            b[(b.index >= '1913-01-01') & (b.index <= '1914-12-01')]).mean()
    p1913 = float(p[p.index.year == 1913].mean())
    us = {}
    for y in YEARS:
        if y <= 1849:
            us[y] = US_WP[y] * 1.0
        elif y <= 1889:
            us[y] = float(a[a.index.year == y].mean() * r_ab) * float(r_bp) / float(p1913) * 100.0
        elif y <= 1912:
            us[y] = float(b[b.index.year == y].mean()) * float(r_bp) / float(p1913) * 100.0
        else:
            us[y] = float(p[p.index.year == y].mean()) / float(p1913) * 100.0
    # Align the 1800-1849 Warren-Pearson block to 1850+ through the 1850 overlap.
    k50 = us[1850] / US_WP[1850]
    for y in range(1800, 1850):
        us[y] = US_WP[y] * k50
    annual['US'] = us
    p13 = p / p1913 * 100.0
    monthly['US'] = p13[p13.index >= '1950-01-01']

    # ===== UK: BOE A47 (2015=100) + ONS GB7S (2015=100) =====
    uk = {}
    wb = __import__('openpyxl').load_workbook(os.path.join(RAW, 'boe_millennium.xlsx'),
                                              read_only=True, data_only=True)
    ws = wb['A47. Wages and prices']
    for r in ws.iter_rows(values_only=True):
        y, v = r[0], r[8]
        if isinstance(y, (int, float)) and isinstance(v, (int, float)):
            uk[int(y)] = float(v)
    wb.close()
    ons = read_ons_gb7s()
    ons_y = ons.groupby(ons.index.year).mean()
    ov = [y for y in range(2008, 2017) if y in uk and y in ons_y.index]
    k = float(np.mean([uk[y] / ons_y[y] for y in ov]))
    for y in range(2008, 2027):
        if y in ons_y.index: uk[y] = float(ons_y[y]) * k
    annual['UK'] = uk
    oecd_uk = read_fred_csv('fred_GBRPROINDMISMEI.csv')   # 2015=100, 1948+
    ov2 = [y for y in range(2008, 2024) if y in ons_y.index]
    k2 = float(np.mean([ons_y[y] / oecd_uk[oecd_uk.index.year == y].mean() for y in ov2]))
    uk_m = pd.concat([oecd_uk * k2, ons])
    uk_m = uk_m[~uk_m.index.duplicated(keep='first')].sort_index()
    monthly['UK'] = uk_m[uk_m.index >= '1950-01-01']

    # ===== DE =====
    de = {}
    ham_a = read_fred_csv('fred_A04054DE00HAMA314NNBR.csv')   # annual Hamburg
    ham_m = read_fred_csv('fred_M04054DE00HAMM314NNBR.csv')   # monthly Hamburg
    for y in range(1851, 1915):
        if y <= 1900 and y in ham_a.index.year.tolist():
            de[y] = float(ham_a[ham_a.index.year == y].mean())
        else:
            sel = ham_m[ham_m.index.year == y]
            if len(sel): de[y] = float(sel.mean())
    base13 = de.get(1913)
    if not base13:
        base13 = float(ham_m[ham_m.index.year == 1913].mean())
    for y in list(de):
        de[y] = de[y] / base13 * 100.0
    de.update(REICHSAMT)
    de.update(DE_1948_57)
    oecd_de = read_fred_csv('fred_DEUPROINDMISMEI.csv')   # 2015=100, 1958+
    k_de = 141.0 / float(oecd_de[oecd_de.index.year == 1958].mean())   # 1958 (1913=100) ~ 141
    for y in range(1958, 2025):
        sel = oecd_de[oecd_de.index.year == y]
        if len(sel): de[y] = float(sel.mean()) * k_de
    dest = read_destatis()   # 2021=100, 2010+
    dest_s = pd.Series(dest)
    dest_s.index = pd.DatetimeIndex([pd.Timestamp(y, m, 1) for y, m in dest_s.index])
    dest_s = dest_s.sort_index()
    ov_d = oecd_de[(oecd_de.index >= '2010-01-01') & (oecd_de.index <= '2024-01-01')]
    ov_d2 = dest_s[(dest_s.index >= '2010-01-01') & (dest_s.index <= '2024-01-01')]
    k_d2 = float((ov_d2 / ov_d.reindex(ov_d2.index).values).mean())
    dest13 = dest_s * k_d2 * k_de   # Destatis -> 1913=100
    for y in range(2010, 2027):
        sel = dest13[dest13.index.year == y]
        if len(sel): de[y] = float(sel.mean())
    annual['DE'] = de
    de_m13 = pd.concat([oecd_de * k_de, dest13])
    de_m13 = de_m13[~de_m13.index.duplicated(keep='first')].sort_index()
    monthly['DE'] = de_m13[de_m13.index >= '1950-01-01']

    # ===== JP =====
    w33 = parse_imes('imes_WPI1933.csv'); w00 = parse_imes('imes_WPI1900.csv'); w87 = parse_imes('imes_WPI1887.csv')
    base3436 = float(np.mean([v for (y, m), v in w33.items() if y in (1934, 1935, 1936)]))
    f33 = 100.0 / base3436
    r31 = float(np.mean([w33[(y, m)] for (y, m) in w33 if y == 1931]) /
                np.mean([w00[(y, m)] for (y, m) in w00 if y == 1931]))
    f00 = f33 * r31
    r01 = float(np.mean([w00[(y, m)] for (y, m) in w00 if y == 1901]) /
                np.mean([w87[(y, m)] for (y, m) in w87 if y == 1901]))
    f87 = f00 * r01
    jp = {}
    for y in YEARS:
        v = [val * f87 for (yy, m), val in w87.items() if yy == y]
        if not v: v = [val * f00 for (yy, m), val in w00.items() if yy == y]
        if not v: v = [val * f33 for (yy, m), val in w33.items() if yy == y]
        if v: jp[y] = float(np.mean(v))
    # IMES (Tokyo WPI) and BOJ national WPI differ in level; splice through a
    # 1945-46 overlap factor (~1.16).
    k_spl = np.mean([JP_BOJ[1945] / jp[1945], JP_BOJ[1946] / jp[1946]]) if 1945 in jp else 1.0
    for y in list(jp):
        if y <= 1944: jp[y] = jp[y] * k_spl
    jp.update(JP_BOJ)
    oecd_jp = read_fred_csv('fred_JPNPROINDMISMEI.csv')   # 2015=100, 1955+
    ov = [y for y in range(1955, 1960) if y in jp]
    k_jp = float(np.mean([jp[y] / oecd_jp[oecd_jp.index.year == y].mean() for y in ov]))
    for y in range(1955, 2025):
        sel = oecd_jp[oecd_jp.index.year == y]
        if len(sel): jp[y] = float(sel.mean()) * k_jp
    boj = parse_boj_tail()
    jp_m = {}
    if boj:
        last = float(oecd_jp[(oecd_jp.index.year == 2024) & (oecd_jp.index.month == 3)].iloc[0]) * k_jp
        for (yy, mm) in sorted(boj):
            last *= 1 + boj[(yy, mm)] / 100.0
            jp_m[(yy, mm)] = last
        for y in range(2024, 2027):
            vals = [v for (yy, mm), v in jp_m.items() if yy == y]
            if vals: jp[y] = float(np.mean(vals))
    annual['JP'] = jp
    jp_m13 = (oecd_jp * k_jp)
    jp_tail = pd.Series({pd.Timestamp(yy, mm, 1): v for (yy, mm), v in jp_m.items()}).sort_index()
    jp_m13 = pd.concat([jp_m13, jp_tail])
    jp_m13 = jp_m13[~jp_m13.index.duplicated(keep='first')].sort_index()
    monthly['JP'] = jp_m13[jp_m13.index >= '1950-01-01']

    # ===== FR =====
    fr = {}
    m = read_fred_csv('fred_M04057FRM360NNBR.csv')
    for y in range(1913, 1940):
        sel = m[m.index.year == y]
        if len(sel): fr[y] = float(sel.mean())
    fr38 = fr[1938]
    # FRB (1938=100, 1940-45) and INSEE (1938=100, 1946+) share the 1946
    # overlap: FRB 1946 level ~648 vs INSEE 606 (both 1938=100).
    k_fb = 606.0 / 648.0
    for y, v in FRB_FR.items():
        fr[y] = v * k_fb / 100.0 * fr38
    for y, v in INSEE_FR.items():
        fr[y] = v / 100.0 * fr38
    oecd_fr = read_fred_csv('fred_FRAPROINDMISMEI.csv')   # 2015=100, 1956+
    insee = read_insee_ippi()   # 2015=100, 2005+
    insee_s = pd.Series(insee)
    insee_s.index = pd.DatetimeIndex([pd.Timestamp(y, m, 1) for y, m in insee_s.index])
    insee_s = insee_s.sort_index()
    ov = oecd_fr[(oecd_fr.index >= '2005-01-01') & (oecd_fr.index <= '2024-01-01')]
    ov2 = insee_s[(insee_s.index >= '2005-01-01') & (insee_s.index <= '2024-01-01')]
    k_fr = float((ov2 / ov.reindex(ov2.index).values).mean())
    fr_m = pd.concat([oecd_fr * k_fr, insee_s])
    fr_m = fr_m[~fr_m.index.duplicated(keep='first')].sort_index()
    # Anchor: French wholesale prices (1938=100) ~3000 in 1956 (1952: 3063,
    # post-war drift ~2%, INSEE long series) -> 1913=100 = 3000/100*653 = 19590.
    fr56 = float(fr_m[fr_m.index.year == 1956].mean())
    k_fr13 = 19590.0 / fr56
    for y in range(1956, 2026):
        sel = fr_m[fr_m.index.year == y]
        if len(sel): fr[y] = float(sel.mean()) * k_fr13
    # 1953-55 gap: linear interpolation on the 1913=100 scale.
    if 1952 in fr and 1956 in fr:
        for y in [1953, 1954, 1955]:
            t = (y - 1952) / 4.0
            fr[y] = fr[1952] * (1 - t) + fr[1956] * t
    annual['FR'] = fr
    monthly['FR'] = fr_m[fr_m.index >= '1950-01-01'] * k_fr13

    return annual, monthly

# ============ 4. Gold-denominated index ============
def gold_index(annual_wpi, fx, gold):
    """gold index = WPI / (FX x gold price), rebased to 1913 = 100."""
    out = {}
    for c, d in annual_wpi.items():
        out[c] = {}
        for y in YEARS:
            v, f, g = d.get(y), fx[c].get(y), gold.get(y)
            if v is None or f is None or g is None or np.isnan(v) or np.isnan(f) or np.isnan(g):
                out[c][y] = np.nan
            else:
                b = (d.get(1913), fx[c].get(1913), gold.get(1913))
                if None in b or np.isnan(b[0]):
                    out[c][y] = np.nan
                else:
                    out[c][y] = (v / (f * g)) / (b[0] / (b[1] * b[2])) * 100.0
    return out

def find_cycles_9ma(s):
    """Peak/trough cycles on the centred 9-yr MA of an annual series:
    local extrema of the MA (>=5 points inside a +-2-yr window), min. 3-yr
    spacing, swing >= 12%, alternating types; a first trough within 5 years of
    the series start is discarded as an edge artefact."""
    m9 = s.rolling(9, center=True, min_periods=9).mean().dropna()
    raw = []
    for y in m9.index:
        window = m9[(m9.index >= y - 2) & (m9.index <= y + 2)]
        if len(window) < 5: continue
        v = float(m9[y])
        is_peak = v >= window.max()
        is_trough = v <= window.min()
        if not (is_peak or is_trough): continue
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
        swing = abs(v - last_v) / max(abs(last_v), 1e-9)
        if swing < 0.12: continue
        pts.append((y, typ, round(v, 1)))
    if not pts and len(m9) > 3:
        pts.append((int(m9.index[0]), 'start', round(float(m9.iloc[0]), 1)))
    return pts

def main():
    gold_a, gold_m = gold_annual_series()
    fx_a, fx_m = fx_annual_series()
    wpi_a, wpi_m = wpi_annual_and_monthly()
    gi = gold_index(wpi_a, fx_a, gold_a)
    countries = ['US', 'UK', 'DE', 'JP', 'FR']

    def rnd(v):
        if v is None:
            return None
        try:
            if np.isnan(v):
                return None
        except TypeError:
            return None
        return round(float(v), 4)

    # ---- Annual table: data/wpi_annual_1800_2026.json ----
    annual = {'gold_usd': {str(y): rnd(gold_a.get(y)) for y in YEARS}}
    for c in countries:
        annual.setdefault('wpi', {})[c] = {str(y): rnd(wpi_a[c].get(y)) for y in YEARS}
        annual.setdefault('fx', {})[c] = {str(y): rnd(fx_a[c].get(y)) for y in YEARS}
        annual.setdefault('goldwpi', {})[c] = {str(y): rnd(gi[c].get(y)) for y in YEARS}
    annual_meta = {
        'title': 'Annual wholesale price index, exchange rate and gold-denominated WPI, 1800-2026',
        'base_year': 1913,
        'goldwpi_units': '1913 = 100 for every country',
        'wpi_units': 'local-currency index on its published base: US/DE/FR 1913=100, '
                     'UK 2015=100, JP 1934-36=100',
        'fx_units': 'local currency per USD (US fixed at 1)',
        'gold_usd_units': 'USD per troy ounce (official price to 1967, then market average)',
        'coverage': {'US': '1800-2026', 'UK': '1800-2026',
                     'DE': '1851-2026 (no data 1945-47)',
                     'JP': '1887-2026', 'FR': '1913-2025'},
    }
    with open(os.path.join(OUT, 'wpi_annual_1800_2026.json'), 'w', encoding='utf-8') as f:
        json.dump({'meta': annual_meta, 'annual': annual}, f, ensure_ascii=False, indent=1)

    # ---- Quarterly 1982+: data/wpi_quarterly_1982_2026.json ----
    qrows = {}
    for c in countries:
        wm, fm, gm = wpi_m[c], fx_m[c], gold_m
        wq = wm[wm.index >= '1982-01-01'].resample('QE').mean()
        fq = fm[fm.index >= '1982-01-01'].resample('QE').mean()
        gq = gm[gm.index >= '1982-01-01'].resample('QE').mean()
        idx = wq.index.intersection(fq.index).intersection(gq.index)
        # Align the quarterly level to the annual 1913=100 scale.
        v = (wq / (fq * gq)).reindex(idx)
        b = (wpi_a[c].get(1913), fx_a[c].get(1913), gold_a.get(1913))
        scale = 100.0 * b[1] * b[2] / b[0] if b[0] else np.nan
        for ts in idx:
            qrows.setdefault(ts, {})[c] = None if np.isnan(v[ts]) else round(float(v[ts] * scale), 3)
    q_list = [{'quarter': f'{ts.year}-Q{(ts.month - 1) // 3 + 1}',
               'US': qrows[ts].get('US'), 'UK': qrows[ts].get('UK'),
               'DE': qrows[ts].get('DE'), 'JP': qrows[ts].get('JP'),
               'FR': qrows[ts].get('FR')} for ts in sorted(qrows)]
    q_meta = {'title': 'Quarterly gold-denominated WPI, 1982-Q1 onwards (monthly PPI / '
                       'FX / gold price averaged per quarter)',
              'units': '1913 = 100',
              'fr_end': '2025-Q4 (INSEE IPPI monthly release paused after 2025-10)'}
    with open(os.path.join(OUT, 'wpi_quarterly_1982_2026.json'), 'w', encoding='utf-8') as f:
        json.dump({'meta': q_meta, 'quarterly': q_list}, f, ensure_ascii=False, indent=1)

    # ---- Cycles on the 9-yr MA of the gold-denominated series ----
    cyc = {}
    for c in countries:
        s = pd.Series({y: gi[c].get(y, np.nan) for y in YEARS}).dropna()
        cyc[c] = [{'year': int(y), 'type': t, 'ma9': v}
                  for (y, t, v) in find_cycles_9ma(s)]
    cyc_meta = {'title': 'Peak/trough cycles of the annual gold-denominated WPI',
                'detection': 'centred 9-yr MA; local extrema with >=5 points in a '
                             '+-2-yr window; min. 3-yr spacing; swing >= 12%; '
                             'alternating types; first trough within 5 years of the '
                             'series start dropped as an edge artefact',
                'caveat': 'the 1922-23 German hyperinflation distorts the 9-yr MA '
                          'around 1918-27; cycle years there carry limited meaning'}
    with open(os.path.join(OUT, 'wpi_cycles_gold.json'), 'w', encoding='utf-8') as f:
        json.dump({'meta': cyc_meta, 'cycles': cyc}, f, ensure_ascii=False, indent=1)
    print('wrote data/wpi_annual_1800_2026.json, data/wpi_quarterly_1982_2026.json, '
          'data/wpi_cycles_gold.json')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Build the gold-denominated WPI series.')
    ap.add_argument('--raw-dir', default=RAW,
                    help='directory with the externally downloaded source files '
                         '(default: <repo>/wpi_raw_data)')
    args = ap.parse_args()
    RAW = os.path.abspath(args.raw_dir)
    main()
