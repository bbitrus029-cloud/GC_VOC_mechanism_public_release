"""Exercise the seven explicit scientific inputs; all fixtures/renders stay external."""
import csv, hashlib, json, os, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipelines'))
from public_io import PublicIO
from figure_scientific_inputs import read, validate

def main():
    if len(sys.argv)!=3: raise ValueError('Usage: test_figure_scientific_inputs.py SEVEN_INPUT_MANIFEST EMPTY_EXTERNAL_OUTPUT')
    root=Path.cwd().resolve();out=Path(sys.argv[2]).resolve();io=PublicIO(sys.argv[1],out)
    spec=json.loads(io.input('workflow/config/figure_scientific_contracts.json').read_text())['inputs']
    original={k:read(io.input(k)) for k in spec};out.mkdir(parents=True)
    checks=[]
    def write(path,rows,fields=None,variant=False):
        with path.open('w',newline='') as f:
            writer=csv.DictWriter(f,fields or list(rows[0]),delimiter='\t',lineterminator='\r\n' if variant else '\n',quoting=csv.QUOTE_ALL if variant else csv.QUOTE_MINIMAL)
            writer.writeheader();writer.writerows(rows)
    def manifest(name,values):
        p=out/(name+'.tsv')
        with p.open('w',newline='') as f:
            writer=csv.DictWriter(f,['logical_input','local_path','sha256'],delimiter='\t',lineterminator='\n');writer.writeheader()
            for k,path in values.items():writer.writerow({'logical_input':k,'local_path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
        return p
    clean={};variant={}
    for k,c in spec.items():
        rows=validate(original[k],c);clean[k]=out/(k+'_clean.tsv');write(clean[k],rows)
        alternate=[]
        for r in reversed(rows):
            v=r.copy()
            for f in c['numeric']:v[f]=format(float(v[f]),'.17e')
            v['binding_file']='IGNORED_OPERATIONAL_METADATA';v['display_x_mm']='999'
            alternate.append(v)
        variant[k]=out/(k+'_equivalent.tsv');write(variant[k],alternate,list(reversed(alternate[0])),True)
        assert validate(read(variant[k]),c)==validate(read(clean[k]),c) or all(
            validate(read(variant[k]),c)[i][f]==rows[i][f] for i in range(len(rows)) for f in c['fields'] if f not in c['numeric'])
        assert hashlib.sha256(clean[k].read_bytes()).hexdigest()!=hashlib.sha256(variant[k].read_bytes()).hexdigest()
        checks.append({'logical_input':k,'case':'clean_and_equivalent_serialization','status':'PASS'})
    cm,vm,empty=manifest('clean_manifest',clean),manifest('equivalent_manifest',variant),manifest('empty_manifest',{})
    def call(label,n,journal,mode,m,success,error=None):
        dest=out/label
        cmd=['workflow/environment/run-r','workflow/pipelines/render_dual_journal_figures.R',str(n),journal,mode,str(m),str(dest)]
        p=subprocess.run(cmd,capture_output=True,text=True);(out/(label+'.log')).write_text(p.stdout+p.stderr)
        if success:assert p.returncode==0 and (dest/'qa/rendered-marks.json').exists(),(label,p.stderr)
        else:
            assert p.returncode!=0 and (not dest.exists() or not any(dest.iterdir())),(label,p.stderr)
            if error:assert error in p.stdout+p.stderr,(label,p.stderr)
        checks.append({'case':label,'command':cmd,'expected_success':success,'returncode':p.returncode,'failed_before_output':not success})
        return dest
    for k,c in spec.items():
        n=int(k.split('_')[1])
        missing=manifest(k+'_missing',{key:p for key,p in clean.items() if key!=k})
        for j in ('JTIM','WJG'):call(k+'_absent_'+j,n,j,'complete',missing,False,'Missing declared input')
        for kind in ('missing_field','wrong_identity','wrong_value','missing_row'):
            rows=[r.copy() for r in original[k]]
            if kind=='missing_field':
                for r in rows:r.pop(c['fields'][0])
            elif kind=='wrong_identity':rows[0][c['keys'][0]]='SYNTHETIC_WRONG_IDENTITY'
            elif kind=='wrong_value':
                field=(c['numeric'] or c['boolean'])[0];rows[0][field]=str(float(rows[0][field])+1) if c['numeric'] else ('FALSE' if rows[0][field]=='TRUE' else 'TRUE')
            else:rows=rows[:-1]
            p=out/(k+'_'+kind+'.tsv');write(p,rows);paths={**clean,k:p};m=manifest(k+'_'+kind+'_manifest',paths)
            for j in ('JTIM','WJG'):call(k+'_'+kind+'_'+j,n,j,'complete',m,False,'Scientific figure input rejected')
    # Actual valid complete routes and equivalent serializations in both journals.
    for n in (2,4,5,8):
        for j in ('JTIM','WJG'):
            a=call(f'complete_{n}_{j}',n,j,'complete',cm,True)
            b=call(f'equivalent_{n}_{j}',n,j,'complete',vm,True)
            for name in ('rendered-marks.json','data-layout.json','panel-layout.json','presentation.json'):
                assert json.loads((a/'qa'/name).read_text())==json.loads((b/'qa'/name).read_text()),(n,j,name)
            assert next(a.glob('*_inspection_192.png')).read_bytes()==next(b.glob('*_inspection_192.png')).read_bytes()
    # Shared loader regression: all earlier bundled-complete / eligible paths.
    for n in range(2,9):
        for j in ('JTIM','WJG'):
            call(f'preserved_{n}_{j}',n,j,'complete' if n in (3,6,7) else 'eligible-panels',empty,True)
    io.guard();(out/'checks.json').write_text(json.dumps({'status':'PASS','checks':checks,'private_fixture_limit':'Existing accepted held analysis values, explicitly supplied; not raw accession reconstruction or redistribution permission.'},indent=2)+'\n')
    print('PASS',len(checks),'scientific input and complete/eligible figure route checks')

def test_preparation():
    if len(sys.argv)!=5:raise ValueError('Usage: test_figure_scientific_inputs.py preparation PROTEIN_MANIFEST COMPOSITION_MANIFEST EMPTY_EXTERNAL_OUTPUT')
    root=Path.cwd().resolve();out=Path(sys.argv[4]).resolve()
    if out.exists() or out==root or root in out.parents:raise ValueError('Use an absent external output')
    out.mkdir(parents=True);spec=json.loads((root/'workflow/config/figure_scientific_contracts.json').read_text());checks=[]
    def write(p,rows,variant=False):
        with p.open('w',newline='') as f:
            fields=list(reversed(rows[0])) if variant else list(rows[0])
            wr=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\r\n' if variant else '\n',quoting=csv.QUOTE_ALL if variant else csv.QUOTE_MINIMAL);wr.writeheader();wr.writerows(rows)
    for kind,source_manifest in zip(('protein','composition'),sys.argv[2:4]):
        original=read(source_manifest)
        for case in ('intact','equivalent','missing','declared_mismatch','content_mismatch','missing_field','shadow','unsafe_output'):
            rows=[r.copy() for r in original];key=rows[0]['logical_input']
            for r in rows:r['local_path']=str((Path(source_manifest).parent/r['local_path']).resolve())
            if case=='missing':rows=rows[1:]
            elif case=='declared_mismatch':rows[0]['sha256']='0'*64
            elif case in ('equivalent','content_mismatch','missing_field'):
                values=read(rows[0]['local_path'])
                if case=='content_mismatch':values[0]['patient_id' if kind=='protein' else 'donor']='SYNTHETIC_WRONG_IDENTITY'
                elif case=='missing_field':
                    field=spec['preparation'][kind]['input_schemas'][key][0]
                    for r in values:r.pop(field)
                p=out/(kind+'_'+case+'_input.tsv');write(p,values,case=='equivalent');rows[0]['local_path']=str(p);rows[0]['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
            elif case=='shadow':rows[0]['logical_input']='workflow/config/figure_scientific_contracts.json'
            m=out/(kind+'_'+case+'_manifest.tsv')
            with m.open('w',newline='') as f:
                wr=csv.DictWriter(f,['logical_input','local_path','sha256'],delimiter='\t',lineterminator='\n');wr.writeheader();wr.writerows(rows)
            dest=root if case=='unsafe_output' else out/(kind+'_'+case)
            cmd=[sys.executable,'workflow/pipelines/prepare_figure_inputs.py',kind,str(m),str(dest)]
            p=subprocess.run(cmd,capture_output=True,text=True);(out/(kind+'_'+case+'.log')).write_text(p.stdout+p.stderr)
            success=case in ('intact','equivalent')
            if success:
                assert p.returncode==0,(kind,case,p.stderr)
                produced=read(dest/'figure_inputs_manifest.tsv');assert len(produced)==(5 if kind=='protein' else 2)
                for r in produced:validate(read(r['local_path']),spec['inputs'][r['logical_input']])
            else:assert p.returncode!=0 and (case=='unsafe_output' or not dest.exists()),(kind,case,p.stderr)
            checks.append({'route':kind,'case':case,'command':cmd,'expected_success':success,'returncode':p.returncode})
    (out/'checks.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n');print('PASS',len(checks),'preparation entry checks')

if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='preparation':test_preparation()
    else:main()
