#!/bin/bash
zcat ${1} | sed 's/^\([0-9XYM]\)/chr\1/g' | bgzip > ${1%.vcf.gz}.chr.vcf.gz
