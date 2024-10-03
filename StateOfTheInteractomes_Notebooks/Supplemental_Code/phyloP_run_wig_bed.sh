#!/bin/bash

for fn in `ls *.wigFix.gz`; do gunzip -c ${fn} | wig2bed - > ${fn}.bed | gzip ; done

