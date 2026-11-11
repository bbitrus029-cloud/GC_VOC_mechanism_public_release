"""Fail-closed comparison regressions; supplied completed runs are read-only."""
import json, os, shutil, subprocess, sys
from pathlib import Path
sys.dont_write_bytecode=True

def main():
    if len(sys.argv)!=5: raise ValueError('Usage: test_second_post_effect_comparison.py B_RUN D_RUN E_RUN EMPTY_EXTERNAL_OUTPUT')
    root=Path.cwd().resolve();out=Path(sys.argv[4]).resolve()
    if root==out or root in out.parents or out.exists(): raise ValueError('Use an absent external test output')
    out.mkdir(parents=True);checks=[]
    expected={
        'B':['results/B_diagnostic.tsv','results/B_donor_values.tsv','diagnostics/B_omissions.tsv','diagnostics/B_bootstrap_draws.tsv','diagnostics/B_g_av.json'],
        'D':['results/D_results.tsv','results/D_patient_deltas.tsv','diagnostics/D_omissions.tsv','diagnostics/D_component_correlations.tsv','diagnostics/D_model_diagnostics.json','diagnostics/D_frozen_reference_moments.tsv','diagnostics/D_influence_summary.tsv'],
        'E':['results/E_results.tsv','results/E_paired_values.tsv','diagnostics/E_bootstrap_draws.tsv','diagnostics/E_omissions.tsv','diagnostics/E_influence_summary.tsv','diagnostics/E_join_population.tsv','diagnostics/E_ratio_identity.tsv','diagnostics/E_aliquot_RNA_bridge.tsv']}
    for stage,original in zip('BDE',sys.argv[1:4]):
        original=Path(original).resolve()
        cases=[('intact',None,None)]
        cases += [('omit_'+str(i),'omit',f) for i,f in enumerate(expected[stage])]
        cases += [('unterminated_'+str(i),'unterminated',f) for i,f in enumerate(expected[stage])]
        cases += [('truncate_'+str(i),'truncate',f) for i,f in enumerate(expected[stage])]
        cases += [(name,name,None) for name in ['empty_scientific','missing_results','missing_diagnostics','wrong_stage_QC','invalid_completion','missing_QC','missing_input_QC','malformed_QC','wrong_stage_file','wrong_stage_mixture','invalid_input_identity']]
        for label,action,file in cases:
            run=out/(stage+'_'+label);shutil.copytree(original,run);report=out/(stage+'_'+label+'_report')
            if action=='omit':(run/file).unlink()
            elif action=='unterminated':(run/file).write_bytes((run/file).read_bytes().rstrip(b'\n'))
            elif action=='truncate':(run/file).write_bytes((run/file).read_bytes()[:max(1,(run/file).stat().st_size//2)])
            elif action=='empty_scientific':
                for d in ['results','diagnostics']:
                    for p in (run/d).iterdir():p.unlink()
            elif action in ('missing_results','missing_diagnostics'):shutil.rmtree(run/action.removeprefix('missing_'))
            elif action in ('wrong_stage_QC','invalid_completion'):
                p=run/'qc'/(stage+'_run.json');s=json.loads(p.read_text());s['stage' if action=='wrong_stage_QC' else 'status']='E' if stage!='E' else 'B';p.write_text(json.dumps(s))
            elif action=='missing_QC':(run/'qc'/(stage+'_run.json')).unlink()
            elif action=='missing_input_QC':(run/'qc/verified_input_identities.tsv').unlink()
            elif action=='malformed_QC':(run/'qc'/(stage+'_run.json')).write_text('{}')
            elif action in ('wrong_stage_file','wrong_stage_mixture'):
                foreign='E' if stage!='E' else 'B';(run/'results'/(foreign+'_results.tsv')).write_text('wrong\tstage\n1\t2\n')
                if action=='wrong_stage_file':(run/expected[stage][0]).unlink()
            elif action=='invalid_input_identity':
                p=run/'qc/verified_input_identities.tsv';lines=p.read_text().splitlines();cells=lines[1].split('\t');cells[-1]='0'*64;lines[1]='\t'.join(cells);p.write_text('\n'.join(lines)+'\n')
            cmd=['workflow/environment/run-r','workflow/pipelines/compare_second_post_effect.R',stage,str(run),str(report)]
            p=subprocess.run(cmd,capture_output=True,text=True);(out/(stage+'_'+label+'.log')).write_text(p.stdout+p.stderr)
            passed=action is None
            if passed:
                assert p.returncode==0,(stage,label,p.stderr)
                receipt=json.loads((report/'comparison.json').read_text())
                assert receipt['status']=='PASS' and receipt['expected_count']==receipt['checked_count']==len(expected[stage])
                assert receipt['expected_outputs']==receipt['checked_outputs']==sorted(expected[stage])
            else:assert p.returncode!=0 and not (report/'comparison.json').exists(),(stage,label,p.stdout,p.stderr)
            checks.append({'stage':stage,'case':label,'expected_pass':passed,'returncode':p.returncode,'command':cmd})
    (out/'checks.json').write_text(json.dumps({'status':'PASS','expected_coverage':expected,'checks':checks},indent=2)+'\n')
    print('PASS',len(checks),'comparison coverage regressions; no calculators executed')

if __name__=='__main__':main()
