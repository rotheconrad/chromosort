#!/usr/bin/env python3
"""Convert standard AGP or RagTag correction AGP to an output/input ledger."""
import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from benchio import write_tsv, read_tsv
FIELDS=['output_id','output_start','output_end','input_id','input_start','input_end','strand']


def read_agp(path, reverse=False):
    rows=[]
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith('#'): continue
        f=line.split('\t')
        if len(f)!=9: raise ValueError('Expected nine AGP columns')
        oid,os,oe=f[0],int(f[1])-1,int(f[2])
        if f[4] in ['N','U']:
            if reverse: raise ValueError('Unexpected gap in correction AGP')
            if int(f[5])!=oe-os: raise ValueError('AGP gap length mismatch')
            rows.append(dict(output_id=oid,output_start=os,output_end=oe,input_id='.',input_start=0,input_end=0,strand='+'))
            continue
        iid,s,e,strand=f[5],int(f[6])-1,int(f[7]),f[8]
        if strand not in ['+','-']: raise ValueError('Unknown AGP orientation cannot receive an invented strand')
        if e-s!=oe-os: raise ValueError('AGP component length mismatch')
        if reverse:
            oid,iid,os,oe,s,e=iid,oid,s,e,os,oe
        rows.append(dict(output_id=oid,output_start=os,output_end=oe,input_id=iid,input_start=s,input_end=e,strand=strand))
    return rows


def compose(rows,parent):
    """Parent maps intermediate records to original input, including gap rows."""
    out=[]
    for r in rows:
        if r['input_id']=='.': out.append(r);continue
        covered=0
        for p in parent:
            if p['output_id']!=r['input_id']: continue
            s,e=max(r['input_start'],int(p['output_start'])),min(r['input_end'],int(p['output_end']))
            if s>=e:continue
            covered+=e-s
            os=r['output_start']+(s-r['input_start'] if r['strand']=='+' else r['input_end']-e)
            if p['input_id']=='.': a=b=0
            elif p['strand']=='+':
                a=int(p['input_start'])+s-int(p['output_start']); b=a+e-s
            else:
                b=int(p['input_end'])-s+int(p['output_start']); a=b-e+s
            out.append(dict(output_id=r['output_id'],output_start=os,output_end=os+e-s,
                            input_id=p['input_id'],input_start=a,input_end=b,
                            strand='+' if p['strand']==r['strand'] else '-'))
        if covered!=r['input_end']-r['input_start']: raise ValueError('Parent ledger does not cover component')
    return sorted(out,key=lambda r:(r['output_id'],r['output_start']))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('agp',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--ragtag-correct',action='store_true');p.add_argument('--parent',type=Path)
    a=p.parse_args();rows=read_agp(a.agp,a.ragtag_correct)
    if a.parent:rows=compose(rows,read_tsv(a.parent))
    write_tsv(a.output,rows,FIELDS)
