"""Read-only current-scope checks used by ./verify; no third-party Python modules."""
import collections, csv, json, math, hashlib
from pathlib import Path
def check(root,table,require,expected):
    read_table=table;cache={}
    def table(path):
        if path not in cache:cache[path]=read_table(path)
        return cache[path]
    R='results/supplementary/second_post_effect/'
    D='results/diagnostics/supplementary/second_post_effect/'
    orig=table('results/integration/original/hypothesis_coverage.tsv')
    formal=[r for r in orig if int(r['frozen_family_size'])>0]
    require(len(formal)==136,'Original 136 formal members plus 10 nonfamily slots')
    require({r['hypothesis_id'] for r in orig}=={r['hypothesis_id'] for r in table('workflow/config/original_hypothesis_estimand_crosswalk.tsv')},'Exact original registered endpoint identities')
    roster=table('workflow/config/second_post_effect_endpoints.tsv');ledger=table(R+'endpoint_ledger.tsv')
    require(len(roster)==10 and [r['hypothesis_id'] for r in roster]==[r['hypothesis_id'] for r in ledger],'Exact PE2 frozen roster order')
    for a,b in zip(roster,ledger):
        require(all(a[k]==b[k] for k in a),'PE2 roster metadata')
        require(b['terminal']=='TRUE' and b['state']!='PENDING','Terminal PE2 ledger')
    for name,n in [('D',6),('E',4)]:
        rows=table(R+name+'_results.tsv');require(len(rows)==n,'Complete PE2 '+name)
        require([r['hypothesis_id'] for r in rows]==[r['hypothesis_id'] for r in roster if r['family']==rows[0]['family']],'Fixed family order '+name)
        finite=[(i,float(r['p_value'])) for i,r in enumerate(rows) if r['p_value']!='NA'];ranked=sorted(finite,key=lambda x:x[1]);q={};current=1.
        for j in range(len(ranked)-1,-1,-1):
            i,p=ranked[j];require(0<=p<=1,'P range');current=min(current,p*n/(j+1));q[i]=current
        for i,r in enumerate(rows):
            require(r['family_size']==str(n) and r['family_closed']=='TRUE','Fixed family size/closure')
            if i in q:
                actual=float(r['q_value']);require(math.isclose(actual,q[i],rel_tol=1e-12,abs_tol=0),'Fixed BH; tiny positive q preserved')
                require((float(r['p_value'])>0)==(actual>0),'No positive P/q rounded to zero')
            else:require(r['q_value']=='NA','Missing q remains missing')
            b=next(x for x in ledger if x['hypothesis_id']==r['hypothesis_id'])
            require(all(b[k]==r[k] for k in ['state','reason','effect','ci_lower','ci_upper','p_value','q_value','biological_n']),'Ledger/result exact values')
    fam=table(R+'family_reconciliation.tsv');require(len(fam)==2,'PE2 two families')
    for r in fam:require(r['frozen_family_size']==r['registered_slots']==r['terminal_slots'] and r['pending_slots']=='0' and r['family_closed']=='TRUE','PE2 family reconciliation')
    disposition=table('publication/RESULT_DISPOSITION.tsv')
    ids=[r['hypothesis_id'] for r in disposition];require(len(ids)==len(set(ids))==198,'All-scope 198 unique slots')
    first=[*table('results/supplementary/external_bulk/all_inferential_results.tsv'),*table('results/supplementary/epithelial/all_results.tsv'),*table('results/integration/supplementary_tissue_protein_hypotheses.tsv')]
    require(set(ids)=={r['hypothesis_id'] for r in orig+first+ledger},'All scopes exact identity union')
    # Preserve the historical 35-row unresolved-source crosswalk verbatim.
    # Current B identity is the separately frozen, already-public GSE150290 roster.
    first_crosswalk=[r for r in table('workflow/config/supplementary_hypothesis_estimand_crosswalk.tsv') if r['module']!='B']
    first_crosswalk += [dict(r,multiplicity_family=r['family'],frozen_family_size=r['family_size']) for r in table('data/processed/GSE150290/inference/complete_hypothesis_roster.tsv')]
    crosswalks=[table('workflow/config/original_hypothesis_estimand_crosswalk.tsv'),first_crosswalk]
    expected_families=table('results/integration/family_reconciliation.tsv')
    for records,cw,label in zip([orig,first],crosswalks,['original','first supplement']):
        require({r['hypothesis_id'] for r in records}=={r['hypothesis_id'] for r in cw},'Exact crosswalk coverage '+label)
        registered={r['hypothesis_id']:r for r in cw};groups=collections.defaultdict(list)
        for r in records:
            c=registered[r['hypothesis_id']];family=r.get('multiplicity_family',r.get('family'));size=r.get('frozen_family_size',r.get('family_size'))
            require(family==c['multiplicity_family'] and size==c['frozen_family_size'],'Frozen membership and family size '+label)
            if int(size)>0:groups[family].append(r)
        require(len(groups)==(29 if label=='original' else 12),'Separate family-instance count '+label)
        for family,members in groups.items():
            require(len(members)==int(registered[members[0]['hypothesis_id']]['frozen_family_size']),'Complete fixed family '+family)
            require(any(r['family']==family and int(r['observed_slots'])==len(members) and r['pending']=='0' for r in expected_families),'Historical family reconciliation '+family)
    catalog=table('workflow/config/second_post_effect_source_objects.tsv')
    groups={r['provenance_group'] for r in catalog}
    for r in table(R+'A_clinical_availability.tsv'):
        require(r['clinical_sources'] in groups and r['population_source'] in groups,'A explicit provenance group')
    for r in catalog:
        objects=json.loads(r['objects_json']);require(bool(objects) and bool(r['access_route']),'A source-object route')
        require(all('/home/' not in o['object'] and '.work/' not in o['object'] for o in objects),'No private A paths')
    backgrounds=table('workflow/config/second_post_effect_background_inputs.tsv')
    require(len({r['logical_input'] for r in backgrounds})==len(backgrounds),'Unique background object identities')

    require('PE2_B_H4_TMM_LOGCPM_Z' not in ids,'B technical diagnostic excluded from inferential slots')
    b=table(R+'B_diagnostic.tsv');require(len(b)==1 and b[0]['family']=='NO_FORMAL_FAMILY' and b[0]['role']=='TECHNICAL_SCORING_SENSITIVITY','B technical identity')
    require([r['donor'] for r in table(R+'B_donor_values.tsv')]==['Pat%02d'%i for i in [*range(1,21),22,23,24]],'B exact 23-pair population/order')
    require(len(table(R+'B_member_values.tsv'))==184 and len(table(D+'B_bootstrap_draws.tsv'))==10000,'B complete measurements/draws')
    require(collections.Counter(r['source'] for r in table(R+'D_patient_deltas.tsv'))=={'GSE122401':80,'TCGA_STAD':33},'D 80/33 patients')
    require(len(table(D+'D_omissions.tsv'))==339 and len(table(D+'E_omissions.tsv'))==320,'Complete influence results')
    join=table(D+'E_join_population.tsv');require(len(join)==80 and sum(int(r['declared_technical_ratios']) for r in join)==130,'E patients versus nested technical repeats')
    e=table(R+'E_paired_values.tsv')
    for g in ['ALDH2','ALDH3A1','ALDH3A2','ALDH7A1']:require([r['patient_id'] for r in e if r['member']==g]==[r['patient_id'] for r in join],'E common ordered population')
    require(len(table(D+'E_bootstrap_draws.tsv'))==10000,'E shared 10000 resamples')
    a=table(R+'A_clinical_availability.tsv');require(len(a)==405,'Complete A aggregate rows')
    for r in a:require(int(r['patient_denominator'])==sum(int(r[k]) for k in ['known_patient_n','unknown_patient_n','unresolved_patient_n']),'A denominator closure')
    for f in ['A_clinical_source_values.tsv','A_population_membership.tsv']:require(not any(r['source'] in {'GSE206785','PDC000614'} for r in table(D+f)),'Held A individual partitions')
    c=table(R+'C_existing_G_and_context_evidence.tsv');require(len(c)==504,'Complete C synthesis')
    for r in c:
        require(r['PE2_new_inferential_slot']=='FALSE','C no new endpoint')
        require(r['source_table'] in expected and hashlib.sha256((root/r['source_table']).read_bytes()).hexdigest()==r['source_sha256'],'C exact public source digest')
        if r['record_type']!='REGISTERED_MISSING_MEMBER_CONTEXT':require(json.loads(r['original_fields_json'])==table(r['source_table'])[int(r['source_row'])-1],'C nested public projection exact')
    k=table(R+'K_IHC_evidence.tsv');require((len(k),len({r['gene'] for r in k}),len({r['antibody'] for r in k}),len({r['source_patient_id'] for r in k}))==(41,4,8,16),'Complete K fixed roster')
    require(all(r['immune_compartment']==r['stromal_compartment']=='NOT_SEPARATELY_ANNOTATED' for r in k),'No invented HPA compartments')
    panels=table('publication/PANEL_MAPPING.tsv');require(len(panels)==32 and len({(r['figure'],r['panel']) for r in panels})==32,'Complete current panels')
    inputs=table('publication/figure_inputs.tsv');require(len({(r['figure'],r['role']) for r in inputs})==len(inputs),'Unique figure roles')
    for r in inputs:
        if r['availability']=='BUNDLED':require(r['canonical_path'] in expected,'Bundled figure identity')
        else:require(not r['canonical_path'] and any(x['logical_input']==r['logical_input'] for x in table('workflow/config/reconstruction_inputs.tsv')),'External figure catalog route')
    scientific=json.loads((root/'workflow/config/figure_scientific_contracts.json').read_text())['inputs']
    external={r['logical_input']:r for r in inputs if r['availability']!='BUNDLED'}
    require(set(scientific)==set(external) and len(external)==7,'Exact seven reconstructed figure contracts')
    reconstruction={r['logical_input']:r for r in table('workflow/config/reconstruction_inputs.tsv')}
    for r in inputs:
        fields=r['schema'].split(';')
        require(not any(k.startswith('binding_') for k in fields),'No private binding fields in reader schemas')
        if r['availability']=='BUNDLED':
            require(fields==list(table(r['canonical_path'])[0]),'Bundled schema describes actual public file')
            require(r['source_sha256']==hashlib.sha256((root/r['canonical_path']).read_bytes()).hexdigest(),'Bundled figure source digest')
        else:
            c=scientific[r['logical_input']];entry=reconstruction[r['logical_input']]
            require(fields==c['fields']==entry['schema'].split(';'),'Consistent reconstructed scientific schema')
            require(entry['binding_rule']==r['binding_rule']=='SCIENTIFIC_CONTENT_CONTRACT' and not entry['authority_sha256'] and not entry['recorded_source_sha256'] and not r['source_sha256'],'Private-view hash is provenance only')
            require(len(c['scientific_sha256'])==64 and c['rows']>0 and set(c['keys'])<=set(fields),'Portable frozen scientific-content authority')
            require(not any(k.startswith('display_') for k in fields),'No operational display coordinates required')
    for p in panels:
        matches=[next(x for x in inputs if x['figure']==p['figure'] and x['role']==role) for role in p['required_roles'].split(';')]
        require((p['bundled_mode']=='COMPLETE_PANEL')==all(x['availability']=='BUNDLED' for x in matches),'Panel completeness cannot overclaim')
    directory=table('publication/TABLE_DIRECTORY.tsv');layout=json.loads((root/'publication/table_layout.json').read_text())
    require(len(directory)==len(layout)==108 and len({r['section_id'] for r in directory})==108,'Current 105 workbook plus 3 reader sections')
    for r,s in zip(directory,layout):require(r['section_id']==s['id'] and r['related_figures']==s['related_figures'] and int(r['eligible_rows'])==len(s['rows']) and int(r['source_rows'])==len(s['rows'])+s['held_rows'],'Section identity/navigation/explicit holds')
    for r in table('publication/table_marks.tsv'):
        require(r['canonical_path'] in expected and table(r['canonical_path'])[int(r['source_row'])-1][r['field']]==r['token'],'Public table source mark')
    allowed={'.md','.tsv','.json','.py','.R','.sh','.yml','.yaml','.txt','.gz'}
    require(all(Path(p).suffix in allowed or p in {'.gitignore','verify','workflow/environment/run-r'} for p in expected),'Default-deny allowed file types')
    dirs={str(p) for f in expected for p in Path(f).parents if str(p)!='.'}
    exceptions={'!/'+p+'/' for p in dirs}|{'!/'+p for p in expected}
    lines=(root/'.gitignore').read_text().splitlines();require('*' in lines and {x for x in lines if x.startswith('!')}==exceptions,'Exact default-deny ignore allowlist')
    return {'second_supplement_hypotheses':10,'second_supplement_families':2,'all_scope_slots':198,'all_scope_family_instances':43,'pooled_FDR':False,'current_panels':32,'table_sections':108}
