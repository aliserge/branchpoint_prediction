#!/usr/bin/env python

import sys


reference_intron = None
reference_intron_data = None

caption_str = ("#rsDI\tchr\tpos\tref\talt\t" +
    'intron\t' +
    'ref_bps\tref_bp_pos\tref_sc_bps\tref_sc_ppt\tref_sc\tref_zsc_bps\tref_zsc_ppt\tref_zsc\t' +
    'alt_bps\talt_bp_pos\talt_sc_bps\talt_sc_ppt\talt_sc\talt_zsc_bps\talt_zsc_ppt\talt_zsc\t' +
    'pos_change\tscore_change\tFactor\tScore\tz_Score\n'
    )
sys.stdout.write(caption_str)
for line in sys.stdin:
    if line.strip() == '' or line[0]=='#' :
        continue
    fields = line.strip().split("\t")
    try:
        intron,variation = fields[0].split(' ')
    except ValueError:
        print(fields)
        exit()
    intron=intron[1:]
    if variation=='ref':
        reference_intron = intron
        reference_intron_data = fields
        continue
    variation = variation.split(':')
    pos_change = int(fields[2]) - int(reference_intron_data[2])
    
    score_ref = float(reference_intron_data[5])
    score_alt = float(fields[5])
    z_score_ref = float(reference_intron_data[-1])
    z_score_alt = float(fields[-1]) 
    
    is_alt = score_alt > score_ref
    factor = 'ALT' if is_alt else 'REF'
    score = score_alt if is_alt else score_ref
    z_score = z_score_alt if is_alt else z_score_ref 
    
    data_fields = (variation + [intron] + reference_intron_data[1:] + fields[1:] + 
        [str(pos_change), str(score_alt - score_ref), factor, str(score), str(z_score)])
    sys.stdout.write('\t'.join(data_fields)+"\n")
    
    
    
