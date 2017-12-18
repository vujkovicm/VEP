#!/bin/bash

grep -w -f ../meta/glist-hg19 > ../output/gene.bed

sed -i "s/ /\t/g" ../output/gene.bed
