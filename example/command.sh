#!/bin/bash

# you must supply a Human Genome reference, such as hg19bis.fa

undr_rover --primer_coords example_coords.txt --primer_sequences example_primers.txt --reference hg19bis.fa --out example.vcf --genotype example_R1.fastq example_R2.fastq 
