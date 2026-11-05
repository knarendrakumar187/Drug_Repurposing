"""Download CTD for local graph reconstruction; KEGG requires your own GMT copy.

No credentials or licences are accepted by this script. Downloads go to a NEW
local file; source data are excluded from version control. CTD terms apply:
https://ctdbase.org/about/legal.jsp
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import urllib.request

URL='https://ctdbase.org/reports/CTD_diseases_pathways.csv.gz'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    if args.out.exists():ap.error('Choose a NEW output file, or use your retained CTD file directly.')
    args.out.parent.mkdir(parents=True,exist_ok=True)
    temp=args.out.with_name(args.out.name+'.partial')
    created=False
    try:
        with urllib.request.urlopen(URL,timeout=120) as response,temp.open('xb') as f:
            created=True
            while block:=response.read(1<<20):f.write(block)
        with gzip.open(temp,'rb') as f:
            for _ in iter(lambda:f.read(1<<20),b''):pass
        if args.out.exists():raise ValueError('Output appeared during download.')
        temp.rename(args.out)
    except BaseException:
        if created:temp.unlink(missing_ok=True)
        raise
    print(json.dumps({'file':args.out.name,'source':URL,'sha256':hashlib.sha256(args.out.read_bytes()).hexdigest()},indent=2))
    print('Keep this retained file; a future CTD update may fail the expected graph-count checks.')


if __name__=='__main__':main()
