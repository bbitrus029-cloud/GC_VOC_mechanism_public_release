#!/usr/bin/env python3
"""Negative public-entry checks using disposable payload copies and explicit inputs."""
import csv, hashlib, json, os, shutil, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipelines'))
from public_io import PublicIO
def main():
    if len(sys.argv) not in (2,3) or (len(sys.argv)==3 and sys.argv[2]!='rendering'):raise ValueError('Usage: test_current_inputs.py EMPTY_EXTERNAL_OUTPUT_DIRECTORY [rendering]')
    root=Path.cwd().resolve();out=Path(sys.argv[1]).resolve()
    io=PublicIO('workflow/config/reconstruction_manifest_template.tsv',out)
    out.mkdir(parents=True,exist_ok=True);fixture=out/'payload';shutil.copytree(root,fixture,ignore=shutil.ignore_patterns('.git'))
    empty=out/'empty.tsv';empty.write_text('logical_input\tlocal_path\tsha256\n');checks=[]
    def call(stage,manifest,output,expected,report=False,render=False,table_render=False):
        command=([sys.executable,'workflow/pipelines/second_post_effect_reporting.py',stage] if report else ['workflow/environment/run-r','workflow/pipelines/second_post_effect.R',stage])+[str(manifest),str(output)]
        if render:command=['workflow/environment/run-r','workflow/pipelines/render_dual_journal_figures.R','6' if stage=='figure6' else '2','JTIM','complete',str(manifest),str(output)]
        if table_render:command=['workflow/environment/run-r','workflow/pipelines/render_data_tables.R','PE2_D','JTIM',str(manifest),str(output)]
        env=os.environ.copy()
        key='GCVOC_RENDER_LIBRARY' if render or table_render else 'GCVOC_PE2_B_LIBRARY' if stage=='B' else 'GCVOC_PE2_DE_LIBRARY'
        if key in env:env['GCVOC_R_LIBRARY']=env[key]
        p=subprocess.run(command,cwd=fixture,env=env,capture_output=True,text=True)
        assert p.returncode!=0 and expected in p.stdout+p.stderr,(stage,p.stdout,p.stderr)
        checks.append({'stage':stage,'expected':expected,'returncode':p.returncode})
    for stage,path,options in [('figure6','publication/figure_inputs.tsv',{'render':True}),('table','publication/table_layout.json',{'table_render':True})]:
        p=fixture/path;data=p.read_bytes()
        p.unlink();call(stage,empty,out/('missing_'+stage),'Bundled input integrity mismatch',**options)
        p.write_bytes(data+b'\n');call(stage,empty,out/('mismatch_'+stage),'Bundled input integrity mismatch',**options);p.write_bytes(data)
        shadow=out/('shadow_'+stage+'.tsv');shadow.write_text('logical_input\tlocal_path\tsha256\n'+path+'\t'+str(p)+'\t'+hashlib.sha256(data).hexdigest()+'\n')
        call(stage,shadow,out/('shadow_'+stage),'cannot be shadowed',**options)
        call(stage,empty,fixture,'outside the immutable release tree',**options)
    if len(sys.argv)==3:
        (out/'checks.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n')
        shutil.rmtree(fixture);print('PASS',len(checks),'rendering entry-point negative checks');return
    for stage,name in [('B','B_member_values.tsv'),('D','D_patient_deltas.tsv'),('E','E_paired_values.tsv')]:
        path='results/supplementary/second_post_effect/'+name;p=fixture/path;data=p.read_bytes()
        p.unlink();call(stage,empty,out/('missing_'+stage),'Bundled input integrity mismatch');p.write_bytes(data+b'\n')
        call(stage,empty,out/('mismatch_'+stage),'Bundled input integrity mismatch');p.write_bytes(data)
        shadow=out/('shadow_'+stage+'.tsv');shadow.write_text('logical_input\tlocal_path\tsha256\n'+path+'\t'+str(p)+'\t'+hashlib.sha256(data).hexdigest()+'\n')
        call(stage,shadow,out/('shadow_'+stage),'cannot be shadowed')
        call(stage,empty,fixture,'outside the immutable release tree')
    for stage,name in [('A','A_clinical_availability.tsv'),('C','C_existing_G_and_context_evidence.tsv'),('K','K_IHC_evidence.tsv')]:
        p=fixture/'results/supplementary/second_post_effect'/name;data=p.read_bytes();p.write_bytes(data+b'\n')
        call(stage,empty,out/('report_mismatch_'+stage),'Input identity mismatch',report=True);p.write_bytes(data)
        call(stage,empty,fixture,'outside the public tree',report=True)
    # Complete figure route must explicitly require its two held protein inputs.
    call('figure2',empty,out/'incomplete_figure','Missing declared input: figure_2_protein_G',render=True)
    # Python catalog shadowing is independently rejected before any reporting.
    path='results/supplementary/second_post_effect/A_clinical_availability.tsv';p=fixture/path
    shadow=out/'shadow_A.tsv';shadow.write_text('logical_input\tlocal_path\tsha256\n'+path+'\t'+str(p)+'\t'+hashlib.sha256(p.read_bytes()).hexdigest()+'\n')
    call('A',shadow,out/'shadow_report','Catalog shadowing',report=True)
    (out/'checks.json').write_text(json.dumps({'status':'PASS','checks':checks},indent=2)+'\n')
    # This entire copy was created by this test and is disposable after checks.
    shutil.rmtree(fixture);print('PASS',len(checks),'public entry-point negative checks')
if __name__=='__main__':main()
