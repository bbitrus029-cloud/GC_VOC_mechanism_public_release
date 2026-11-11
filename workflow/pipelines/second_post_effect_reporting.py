#!/usr/bin/env python3
"""A/C/K reporting from explicit public evidence; no inference or network."""
import collections, json, sys, re
sys.dont_write_bytecode=True
from public_io import PublicIO
R='results/supplementary/second_post_effect/'
D='results/diagnostics/supplementary/second_post_effect/'
def c_display_construct(source):
    missing=(None,'','NA','NA_NOT_APPLICABLE','NONE')
    if source.get('member') not in missing:return source['member']
    member=source.get('members','')
    if str(source.get('member_count'))=='1' and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',member):return member
    return next((source[k] for k in ['construct','construct_id','program','scientific_program'] if source.get(k) not in missing),'SEE_ORIGINAL_FIELDS')
def main():
    if len(sys.argv)!=4 or sys.argv[1] not in {'A','C','K'}:raise ValueError('Usage: second_post_effect_reporting.py A|C|K INPUT_MANIFEST.tsv OUTPUT_DIRECTORY')
    stage=sys.argv[1];io=PublicIO(sys.argv[2],sys.argv[3])
    if stage=='A':
        a=io.table(R+'A_clinical_availability.tsv');membership=io.table(D+'A_population_membership.tsv');values=io.table(D+'A_clinical_source_values.tsv')
        by=collections.defaultdict(list)
        for r in values:by[(r['source'],r['patient_id'],r['field'])].append(r['raw_value'])
        unknown={'','na','n/a','unknown','not reported','not available','nan'};count=0
        for r in a:
            if r['source'] in {'PDC000614','GSE206785'}:continue
            pp=[p for p in membership if all(p[k]==r[k] for k in ['source','population','inclusion'])]
            known=missing=conflict=0;labels=collections.Counter()
            for p in pp:
                v={s for s in by[(p['source'],p['patient_id'],r['field'])] if s.strip().lower() not in unknown}
                if len(v)>1:conflict+=1;labels['UNRESOLVED_MULTIPLE_SPECIMEN_VALUES']+=1
                elif v:known+=1;labels[next(iter(v))]+=1
                else:missing+=1;labels['UNKNOWN_OR_NOT_REPORTED']+=1
            calculated=[len(pp),known,missing,conflict]
            assert calculated==[int(r[k]) for k in ['patient_denominator','known_patient_n','unknown_patient_n','unresolved_patient_n']]
            assert dict(labels)==json.loads(r['source_label_counts']);count+=1
        io.write('A_clinical_availability.tsv',a)
        report={'eligible_aggregate_rows_rebuilt':count,'held_aggregate_rows_reused':len(a)-count,'detailed_held_sources':['PDC000614','GSE206785'],'inference':False}
    elif stage=='C':
        c=io.table(R+'C_existing_G_and_context_evidence.tsv');assert len(c)==504
        for r in c:
            assert r['PE2_new_inferential_slot']=='FALSE'
            if r['record_type']=='REGISTERED_MISSING_MEMBER_CONTEXT':continue
            source=io.table(r['source_table']);original=source[int(r['source_row'])-1]
            assert original==json.loads(r['original_fields_json'])
            if r['record_type']=='REUSED_ACCEPTED_RESULT':assert r['construct']==c_display_construct(original)
        single=[r for r in c if r['source']=='GSE122401' and r['record_type']=='REUSED_ACCEPTED_RESULT' and json.loads(r['original_fields_json']).get('member')]
        # The existing corrected singleton label is the scientific display identity.
        for r in single:
            d=json.loads(r['original_fields_json'])
            if d.get('member') in {'ACSL4','LPCAT3','ALOX5','ALOX12','ALOX15','ALOX15B','POR','ALDH2','ALDH3A1','ALDH3A2','ALDH7A1'}:assert r['construct']==d['member']
        io.write('C_existing_G_and_context_evidence.tsv',c);io.write('C_background_dictionary.tsv',io.table(R+'C_background_dictionary.tsv'))
        report={'rows':504,'reused_results':454,'detection_coverage':49,'missing_member_context':1,'new_endpoints':0,'singleton_correction':'2026-11-08 retained'}
    else:
        k=io.table(R+'K_IHC_evidence.tsv');assert len(k)==41 and len({r['gene'] for r in k})==4 and len({r['antibody'] for r in k})==8 and len({r['source_patient_id'] for r in k})==16
        assert all(r['immune_compartment']=='NOT_SEPARATELY_ANNOTATED' and r['stromal_compartment']=='NOT_SEPARATELY_ANNOTATED' for r in k)
        roster=io.table('provenance/registration/second_post_effect_supplement_2026-11-07/HPA_IHC_FIXED_SOURCE_ROSTER.tsv')
        assert len(roster)==41 and {r['image_id'] for r in k}=={r['image_id'] for r in roster}
        io.write('K_IHC_evidence.tsv',k)
        report={'release':'HPA 25.1','records':41,'genes':4,'antibodies':8,'distinct_source_patient_ids':16,'images_redistributed':False,'RNA_analysis':False,'inference':False,'qualification':'Categorical observations and antibody disagreement retained; RNA-informed validation is a source qualification, not independent replication.'}
    io.guard();io.write('report.json',{'status':'PASS','stage':stage,**report});print(json.dumps(report))
if __name__=='__main__':main()
