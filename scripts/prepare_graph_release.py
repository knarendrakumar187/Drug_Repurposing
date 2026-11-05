"""Verify and extract the published Paper 1 ZIP for a separate Paper 2 import.

Does not connect to Neo4j or overwrite a database. Optionally print the native
import command for the extracted CSVs. The ZIP digest pins the starting release.
"""
import argparse
import csv
import hashlib
import io
from pathlib import Path
import shutil
import tempfile
import zipfile

EXPECTED_SHA256='62e67bcccce63520262c556ef5bef19a4155e5e2db744fb1f51bdfd44f45d2d6'
RELEASE_URL='https://github.com/nissankararaoThanvitha/liver_graph/releases/tag/graph-export-2026-11-05'


def extract(archive,out):
    h=hashlib.sha256()
    with archive.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    if h.hexdigest()!=EXPECTED_SHA256:raise ValueError('ZIP does not match the verified Paper 1 release: '+RELEASE_URL)
    if out.exists():raise ValueError('Choose a NEW extraction directory.')
    out.parent.mkdir(parents=True,exist_ok=True)
    temp=Path(tempfile.mkdtemp(prefix=out.name+'.partial-',dir=out.parent))
    try:
        with zipfile.ZipFile(archive) as z:
            for entry in z.infolist():
                if Path(entry.filename).is_absolute() or '..' in Path(entry.filename).parts:raise ValueError('Unsafe ZIP path.')
            rows=list(csv.DictReader(io.StringIO(z.read('MANIFEST.csv').decode())))
            for row in rows:
                h=hashlib.sha256()
                with z.open(row['file']) as f:
                    for block in iter(lambda:f.read(1<<20),b''):h.update(block)
                if h.hexdigest()!=row['sha256']:raise ValueError('Data checksum mismatch: '+row['file'])
            z.extractall(temp)
        if out.exists():raise ValueError('Output appeared during extraction.')
        temp.rename(out)
        return len(rows)
    except BaseException:
        shutil.rmtree(temp)
        raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--archive',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();print('Verified and extracted',extract(args.archive,args.out),'data files into',args.out)
    print('Import into a NEW, separate Paper 2 database using the enclosed IMPORT.md. Do not overwrite Paper 1.')


if __name__=='__main__':main()
