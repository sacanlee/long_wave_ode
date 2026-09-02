# -*- coding: utf-8 -*-
"""
Download the Destatis producer-price statistical report (the xlsx that backs
the tail of the German WPI series, table 61241-b01, 2021=100, 2010-2026).

The Destatis server is very slow (~50 KB/s, no content-length header), so
plain curl downloads get truncated. This script streams with python requests,
retries up to 4 times and verifies the xlsx (ZIP) integrity.

Usage:
    python scripts/wpi_dl_destatis.py [--raw-dir <dir>]

Writes <raw-dir>/destatis_ppi3.xlsx (default: <repo>/wpi_raw_data/).
"""
import os, sys, argparse, zipfile, time
import requests

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
URL = ('https://www.destatis.de/DE/Themen/Wirtschaft/Preise/'
       'Erzeugerpreisindex-gewerbliche-Produkte/Publikationen/'
       'Downloads-Erzeugerpreise/'
       'statistischer-bericht-erzeugerpreise-2170200261075.xlsx'
       '?__blob=publicationFile&v=3')

def main(raw_dir):
    os.makedirs(raw_dir, exist_ok=True)
    out = os.path.join(raw_dir, 'destatis_ppi3.xlsx')
    for attempt in range(4):
        try:
            r = requests.get(URL, timeout=(30, 900), stream=True)
            total = 0
            with open(out, 'wb') as f:
                for chunk in r.iter_content(1 << 16):
                    f.write(chunk)
                    total += len(chunk)
            print('attempt', attempt, 'downloaded', total, 'bytes')
            try:
                z = zipfile.ZipFile(out)
                print('ZIP OK', len(z.namelist()), 'entries')
                return
            except Exception as e:
                print('zip integrity fail:', e)
                time.sleep(3)
        except Exception as e:
            print('attempt', attempt, 'err:', e)
            time.sleep(3)
    sys.exit('download failed after 4 attempts')

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Download the Destatis PPI statistical report xlsx.')
    ap.add_argument('--raw-dir', default=os.path.join(REPO, 'wpi_raw_data'),
                    help='target directory for the downloaded file '
                         '(default: <repo>/wpi_raw_data)')
    args = ap.parse_args()
    main(os.path.abspath(args.raw_dir))
