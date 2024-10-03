#!/bin/bash

mkdir -p perGene

for name in `cut -f4 genes.v46.basic.bed | sort | uniq`; do grep -F ${name} genes.v46.basic.bed > perGene/${name}.bed; done
