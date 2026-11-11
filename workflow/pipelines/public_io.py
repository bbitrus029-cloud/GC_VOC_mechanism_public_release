"""Explicit, immutable public inputs and safe external outputs (standard library)."""
import csv, hashlib, json, os
from pathlib import Path
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
def rows(path):
    with Path(path).open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
class PublicIO:
    def __init__(self,manifest,output):
        self.root=Path.cwd().resolve();p=self.root/'PUBLIC_MANIFEST.tsv'
        if not p.is_file():raise ValueError('Run from the public release root')
        self.pub={r['path']:r for r in rows(p)};self.checked={p:sha(p)}
        self.alias={r['logical_input']:r['bundled_path'] for r in rows(self.input('workflow/config/bundled_input_aliases.tsv'))}
        self.catalog=rows(self.input('workflow/config/reconstruction_inputs.tsv'))
        raw=Path(output).expanduser().absolute()
        if any(x.is_symlink() for x in (raw,*raw.parents)):raise ValueError('Unsafe output: symbolic link')
        self.out=raw.resolve()
        if self.out==self.root or self.root in self.out.parents:raise ValueError('Output must be outside the public tree')
        if self.out.exists() and (not self.out.is_dir() or any(self.out.iterdir())):raise ValueError('Output must be absent or empty')
        self.mapping={};manifest=Path(manifest).resolve(strict=True);mr=rows(manifest)
        with manifest.open() as f:
            if next(csv.reader(f,delimiter='\t'))!=['logical_input','local_path','sha256']:raise ValueError('Invalid input manifest schema')
        self.checked[manifest]=sha(manifest)
        for r in mr:
            key=r['logical_input'];q=Path(key)
            if q.is_absolute() or '..' in q.parts or key in self.pub or key in self.alias or key in self.mapping:raise ValueError('Catalog shadowing or invalid logical input')
            p=Path(r['local_path']).expanduser();p=(manifest.parent/p).resolve(strict=True)
            if self.out in p.parents:raise ValueError('Input cannot be in output')
            if sha(p)!=r['sha256']:raise ValueError('Declared input identity mismatch')
            self.mapping[key]=(p,r['sha256'])
    def input(self,key):
        key=getattr(self,'alias',{}).get(key,key)
        if key in self.pub:
            p=self.root/key;expected=self.pub[key]['sha256']
            if p.is_symlink() or p.resolve()!=p or p.stat().st_nlink!=1:raise ValueError('Unsafe bundled input')
        elif key in getattr(self,'mapping',{}):
            p,expected=self.mapping[key]
            for r in self.catalog:
                if r['logical_input']==key and r['binding_rule']=='EXACT_ACQUISITION' and r['authority_sha256']!=expected:raise ValueError('Authoritative input identity mismatch')
        else:raise ValueError('Missing explicit input: '+key+'; no fallback or download')
        if sha(p)!=expected:raise ValueError('Input identity mismatch: '+key)
        self.checked[p]=expected;return p
    def table(self,key):return rows(self.input(key))
    def guard(self):
        for p,h in self.checked.items():
            if sha(p)!=h:raise ValueError('Input changed during execution')
    def write(self,name,values):
        self.guard();p=self.out/name
        if Path(name).is_absolute() or '..' in Path(name).parts or p.exists():raise ValueError('Unsafe or existing output')
        p.parent.mkdir(parents=True,exist_ok=True);stage=p.with_name(p.name+'.stage')
        if stage.exists() or stage.is_symlink():raise ValueError('Unresolved staged output')
        if name.endswith('.json'):
            stage.write_text(json.dumps(values,indent=2,ensure_ascii=False)+'\n');json.loads(stage.read_text())
        else:
            with stage.open('w',newline='') as f:
                w=csv.DictWriter(f,list(values[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(values)
            if rows(stage)!=[{k:str(v) for k,v in r.items()} for r in values]:raise ValueError('Persisted TSV mismatch')
        os.replace(stage,p)
