#!/bin/bash
zcat ${1} | sed 's/>\([0-9XYMT]\+\)/>chr\1/g' | bgzip > ${1%.fa.gz}.chr.fa.gz
