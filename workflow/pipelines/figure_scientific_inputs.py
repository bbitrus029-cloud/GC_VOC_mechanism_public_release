"""Portable scientific projections; no private bindings or historical byte hashes.

The frozen digest commits to typed scientific values and identities, not TSV
formatting. Finite numbers use Python float.hex (IEEE-754 binary64); booleans
accept TRUE/FALSE; identifiers are exact UTF-8. Rows sort by declared identity.
No rounding, fitting, source search, download or participant data is embedded.
"""
import csv, hashlib, json, math, sys
from pathlib import Path
sys.dont_write_bytecode = True

def read(path):
    with Path(path).open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ValueError('Invalid scientific input header')
        rows = list(reader)
    if not rows or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError('Empty or truncated scientific input')
    return rows

def project(rows, contract):
    fields = contract['fields']; numeric = contract['numeric']; boolean = contract['boolean']
    if any(not set(fields) <= set(r) for r in rows):
        raise ValueError('Missing required scientific fields')
    result = [{k:r[k] for k in fields} for r in rows]
    result.sort(key=lambda r: tuple(r[k] for k in contract['keys']))
    keys = [tuple(r[k] for k in contract['keys']) for r in result]
    if len(keys) != len(set(keys)) or any(any(not x or x == 'NA' for x in key) for key in keys):
        raise ValueError('Duplicate or missing scientific identity')
    typed = []
    for r in result:
        row = []
        for k in fields:
            v = r[k]
            if k in numeric:
                x = float(v)
                if not math.isfinite(x): raise ValueError('Nonfinite scientific value: '+k)
                row.append(['number', (0.0 if x == 0 else x).hex()])
            elif k in boolean:
                if str(v).upper() not in ('TRUE','FALSE'): raise ValueError('Invalid qualification flag')
                row.append(['boolean', str(v).upper()])
            else:
                if not v or v == 'NA': raise ValueError('Missing scientific identity/scale: '+k)
                row.append(['text', v])
        typed.append(row)
    canonical = json.dumps([fields,typed],ensure_ascii=False,separators=(',',':')).encode('utf-8')
    return result, hashlib.sha256(canonical).hexdigest()

def validate(rows, contract):
    result, digest = project(rows, contract)
    if len(result) != contract['rows'] or digest != contract['scientific_sha256']:
        raise ValueError('Frozen scientific identity/content mismatch (population, pairing, members, units or values)')
    return result

def prepare(kind, tables):
    """Project existing public analysis outputs; never rerun their models."""
    H = ['ALDH2','ALDH3A1','ALDH3A2','ALDH7A1']
    if kind == 'protein':
        pairs, members, bridge = (tables[k] for k in ('figure_source_protein_pairs','figure_source_protein_members','figure_source_bridge_pairs'))
        out = {}
        for g in ('G','H'):
            out['figure_2_protein_'+g] = [r for r in pairs if r['endpoint_id']=='PDC000614_protein_'+g]
            out['figure_5_'+g+'_pairs'] = [r for r in bridge if r['endpoint_id']=='CPTAC_bridge_'+g]
        selected = [r for r in members if r['population_reference_id']=='protein_H_program' and r['member'] in H and r['layer']=='protein']
        groups = {}
        for r in selected:
            key=(r['member'],r['patient_id']); tissues=groups.setdefault(key,{})
            if r['tissue'] not in ('NAT','Tumor') or r['tissue'] in tissues: raise ValueError('Invalid member pairing')
            tissues[r['tissue']]=r
        values=[]
        for (member,patient),tissues in sorted(groups.items()):
            if set(tissues)!= {'NAT','Tumor'}: raise ValueError('Incomplete member pair')
            n,t=tissues['NAT'],tissues['Tumor']
            if n['input_scale'] != t['input_scale']: raise ValueError('Inconsistent paired units')
            values.append(dict(population_reference_id='protein_H_program',member=member,patient_id=patient,
                NAT_value=n['value'],Tumor_value=t['value'],delta=format(float(t['value'])-float(n['value']),'.17g'),
                NAT_sample_id=n['sample_id'],Tumor_sample_id=t['sample_id'],input_scale=n['input_scale']))
        out['figure_4_protein_changes']=values
        return out
    if kind == 'composition':
        rows=tables['figure_source_composition']; donors={}
        for r in rows:
            d={k:r[k] for k in ('donor','primary_100','strict_200','count_aware')}
            if r['donor'] in donors and donors[r['donor']] != d: raise ValueError('Inconsistent donor eligibility')
            donors[r['donor']]=d
        return {'figure_8_composition':rows,'figure_8_donors':list(donors.values())}
    raise ValueError('Unknown preparation route')

if __name__ == '__main__':
    if len(sys.argv)!=4: raise SystemExit('Usage: figure_scientific_inputs.py CONTRACT.json LOGICAL_INPUT INPUT.tsv')
    contract=json.loads(Path(sys.argv[1]).read_text())['inputs'][sys.argv[2]]
    values=validate(read(sys.argv[3]),contract)
    print(json.dumps({'status':'PASS','logical_input':sys.argv[2],'rows':len(values),'scientific_sha256':contract['scientific_sha256']}))
