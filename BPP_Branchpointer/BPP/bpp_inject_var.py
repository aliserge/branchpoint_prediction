#!/usr/bin/env python3

import sys
from collections import namedtuple, defaultdict
from Bio import Seq,SeqIO

IntronHeader = namedtuple('IntronHeader', ['intron_id','chrom','start','stop','direction'])

Variation = namedtuple('Variation', ['rsid','chrom','pos','ref','alt'])

# Read infromation about intron fro fasta header
def parse_intron_header(header):
    intron_id,_,chrom,pos = header.split(':')
    direction = pos[-2]
    start,stop = [ int(i) for i in pos[:-3].split('-') ]
    return IntronHeader(intron_id,chrom,start,stop,direction)

# load dictionary of intrones, collected by chromosome 
def load_intrones(introns_file):
    introns_dict = defaultdict(list)
    for intron in SeqIO.parse(introns_file, "fasta"):
        header = parse_intron_header(intron.id)
        seq_len = header.stop-header.start
        if seq_len != len(intron.seq):
            raise ValueError(f"sequence length {seq_len} don't match coordinates form intron {header}\n")
        introns_dict[header.chrom].append((header,intron.seq))
    return introns_dict
    
def parse_varition_str(variation_str):
    rsid,chrom,pos,_,ref,alt = variation_str.strip().split('\t')
    return Variation(rsid,chrom,int(pos),ref,alt )
    
def get_next_variation(variations_file):
    variation_str = variations_file.readline()
    if not variation_str:
        return None
    return parse_varition_str(variation_str)


def find_intron(introns_list, chrom, pos, start_position=0):
    i = start_position
    while i<len(introns_list):
        header, seq = introns_list[i]
        if header.start < pos <= header.stop:
            return header, seq, i
        if header.start > pos:
            break
        i+=1
    return None, None, i

intrones_file_name=sys.argv[1]
variations_file_name=sys.argv[2]


introns_file = open(intrones_file_name)
introns_dict = load_intrones(introns_file)
sys.stderr.write("Introns loaded\n")

parsed_intrones = set()
variations_file = open(variations_file_name)

var_counter = 0
filtered_variations = 0
non_intronic_variations = 0
error_variations = 0
processed_variations = 0

current_intron = 0
current_chrom = ''

while True:
    variation = get_next_variation(variations_file)
    if variation is None:
        break
        
    var_counter += 1
    if var_counter % 100000 == 0:
         sys.stderr.write(f"LOG: processed {var_counter} variations\n")
        
    if not variation.chrom in introns_dict:
        sys.stderr.write("Variation chromosome not found in intrones list\n")
        sys.stderr.write(str(variation)+'\n')
        error_variations += 1
        continue
    if variation.chrom != current_chrom:
        current_intron = 0
        current_chrom = variation.chrom
        
    header, seq, current_intron = find_intron(introns_dict[variation.chrom], variation.chrom, variation.pos, current_intron)
    if header is None:
        non_intronic_variations +=1
        continue
    
    bed_pos = variation.pos -1 # converting to BED notation
    diff_pos = bed_pos - header.start
    
    mod_seq = seq
    if header.direction == '-':
        mod_seq = mod_seq.reverse_complement()    
    if (diff_pos >= len(mod_seq)) or (diff_pos < 0):
        raise RuntimeError(f"Varition position {diff_pos} not matches intron length {header}\n{header}\n{variation}\n")
    
    ref_letter = str(mod_seq[diff_pos:diff_pos+len(variation.ref)])
    ref_letter_len = len(ref_letter)
    if ref_letter_len != len(variation.ref):
        if ref_letter in variation.ref:
            sys.stderr.write(f"Variation REF truncated by the end of intron {ref_letter}\n{header}\n{variation}\n")
        else:
            raise RuntimeError(f"Varition refrence sequence {ref_letter} lenght {len(ref_letter)} not matches intron length {len(variation.ref)}\n{header}\n{variation}\n")
    if ref_letter != variation.ref:
        sys.stderr.write(f"Variation REF not matches reference letter {ref_letter}\n{header}\n{variation}\n")
        filtered_variations += 1
        continue
    
    if not header in parsed_intrones:
        parsed_intrones.add(header)
        record = SeqIO.SeqRecord(seq, id=':'.join([str(i) for i in header]),
            description='ref')
        SeqIO.write(record, sys.stdout, "fasta")
    
    mod_seq = list(mod_seq)
    mod_seq[diff_pos:diff_pos+ref_letter_len] = variation.alt
    mod_seq = Seq.Seq(''.join(mod_seq).replace('N','G'))
    
    if header.direction == '-':
        mod_seq = mod_seq.reverse_complement()  
    
    record = SeqIO.SeqRecord(mod_seq, id=':'.join([str(i) for i in header]),
        description=':'.join([str(i) for i in variation]))
    SeqIO.write(record, sys.stdout, "fasta")
    processed_variations += 1
    
sys.stderr.write(f'Recieved {var_counter} variations\nProcessed {processed_variations}, {non_intronic_variations} non-intronic, {filtered_variations} non-matching reference, {error_variations} errors\n')

