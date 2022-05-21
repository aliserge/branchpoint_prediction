#!/usr/bin/env python3

import sys
from collections import namedtuple

min_intron_length = 50

previous = None

GtfRecord = namedtuple('GtfRecord', ['seqid', 'source',	'type', 'start', 'stop', 'score' ,	'strand', 'phase' , 'attributes'])

def parse_attributes(attr_line):
	attr_dict={}
	for attr in attr_line.split(';'):
		attr=attr.strip()
		if attr=='':
			continue
		#sys.stderr.write('attr'+attr+'\n')
		key, value = attr.strip().split(' ')
		attr_dict[key]=value.strip('"')
	return attr_dict

line_number=0
intron_number=0
for line in sys.stdin:
	line_number +=1
	if line.startswith('#'):
		sys.stdout.write(line)
		continue
	current = GtfRecord(*line.strip().split('\t'))
	if current.type != 'exon':
		previous = None
		intron_number=0
		continue
	
	if previous is None:
		previous = current
		intron_number=1
		continue
		
	if current.type == previous.type:
		start = int(previous.stop)+1
		stop  = int(current.start)-1
		if stop - start >= min_intron_length:
			attr=parse_attributes(current.attributes)
			intron = GtfRecord(current.seqid, current.source,	f'{attr["transcript_id"]}.intron_{intron_number}', 
				str(start), str(stop), current.score, current.strand, current.phase, current.attributes)
			sys.stdout.write('\t'.join(intron)+'\n')
			intron_number += 1
	previous = current
