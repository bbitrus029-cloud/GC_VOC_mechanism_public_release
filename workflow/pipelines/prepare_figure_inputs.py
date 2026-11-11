"""Prepare held scientific views from explicit existing analysis-result tables."""
import json, sys
from pathlib import Path
sys.dont_write_bytecode=True
from public_io import PublicIO, sha
from figure_scientific_inputs import prepare, validate

def main():
    if len(sys.argv)!=4 or sys.argv[1] not in ('protein','composition'):
        raise ValueError('Usage: prepare_figure_inputs.py protein|composition INPUT_MANIFEST.tsv EMPTY_EXTERNAL_OUTPUT_DIRECTORY')
    kind=sys.argv[1];io=PublicIO(sys.argv[2],sys.argv[3])
    spec=json.loads(io.input('workflow/config/figure_scientific_contracts.json').read_text())
    tables={key:io.table(key) for key in spec['preparation'][kind]['inputs']}
    for key,fields in spec['preparation'][kind]['input_schemas'].items():
        if not tables[key] or any(not set(fields)<=set(row) for row in tables[key]):
            raise ValueError('Missing analysis-output scientific fields: '+key)
    projected=prepare(kind,tables)
    # Complete scientific checks precede any output, including partial batches.
    checked={key:validate(value,spec['inputs'][key]) for key,value in projected.items()}
    for key,value in checked.items():io.write(key+'.tsv',value)
    manifest=[{'logical_input':key,'local_path':str(io.out/(key+'.tsv')),'sha256':sha(io.out/(key+'.tsv'))} for key in checked]
    io.write('figure_inputs_manifest.tsv',manifest)
    io.write('preparation.json',{'status':'PASS','route':kind,'scientific_rows':{k:len(v) for k,v in checked.items()},'new_inference':False,'historical_private_byte_hash_required':False})
    io.guard();print('PASS',kind,'scientific projections:',len(checked))

if __name__=='__main__':main()
