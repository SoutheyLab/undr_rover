# UNDR ROVER - Unmapped primer directed read overlap variant caller

Version: 2.0.0

## Authors:

 * Bernard J Pope (2,3)
 * Tú Nguyen-Dumont (1)
 * Daniel J Park (1)
 * Roger Li (2)
 * Edmund Lau (2)
 * Peter Georgeson (2)

1. Genetic Epidemiology Laboratory, Department of Pathology, The University of Melbourne.
2. Victorian Life Sciences Computation Initiative (VLSCI).
3. Department of Computing and Information Systems, The University of Melbourne.
         
* [Web page](https://github.com/bjpop/undr_rover)
* [User Documentation](http://bjpop.github.io/undr_rover/)


## License

[The BSD 3-Clause License](http://opensource.org/licenses/BSD-3-Clause)


## Requirements

  * Python 2.7
  * [PyVCF](https://pypi.python.org/pypi/PyVCF)
  * [Pyfaidx](https://pypi.python.org/pypi/pyfaidx)
  * [Biopython](https://pypi.python.org/pypi/biopython)
  * [SciPy](http://www.scipy.org/)

See the [User Documentation](http://bjpop.github.io/undr_rover/) for information
about installation.

## General Description

Undr Rover enables the user to call variants from decompressed FASTQ files 
directly without having to go through a mapping step. 

Reads are organised into blocks, which are initialised from information about
known primer sequences. For each block, we record the primer's DNA sequence, 
coordinates and reference insert sequence. Once the reads have been assigned
to blocks (by matching the bases in their primer region to a known primer), 
we call variants for each block one at a time. 

The approach Undr Rover takes when processing reads to call variants is to 
initially assume that all reads contain only single nucleotide variants. If 
we detect two single nucleotide variants in a row during this variant calling, 
Undr Rover immediatly ceases this step and instead does a full gapped alignment 
of the read against the reference insert sequence to detect possible insertions 
and deletions. 

Only variants detected in both reads of a read-pair are considered relevant. 
User-defined thresholds determine the minimum number and proportion of 
read-pairs a variant must be observed to be considered a 'PASS'. Coverage is 
also reported for each block to allow for regions which require additional 
screening to be easily identified.

## Command Line Usage

See the [User Documentation](http://bjpop.github.io/undr_rover/) for more information
about how to use UNDR ROVER including a worked example.

```
usage: undr_rover [-h] [--version] --primer_coords PRIMER_COORDS
                  --primer_sequences FILE
                  [--primer_prefix_size PRIMER_PREFIX_SIZE]
                  [--kmer_length KMER_LENGTH]
                  [--kmer_threshold KMER_THRESHOLD]
                  [--primer_bases PRIMER_BASES] [--proportionthresh N]
                  [--absthresh N] [--qualthresh N] [--overlap OVERLAP]
                  [--max_variants N] --reference FILE [--id_info FILE] --out
                  FILE [--log FILE] [--coverdir COVERDIR] [--fast]
                  [--snvthresh N] [--genotype] [--ploidy {1,2}]
                  [--error ERROR]
                  fastqs [fastqs ...]

Find variants from fastqs via a mapping-free approach.

positional arguments:
  fastqs                Fastq files containing reads.

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --primer_coords PRIMER_COORDS
                        Primer coordinates in TSV format.
  --primer_sequences FILE
                        Primer base sequences as determined by a primer
                        generating program.
  --primer_prefix_size PRIMER_PREFIX_SIZE
                        Size of primer key for blocks in dictionary.
  --kmer_length KMER_LENGTH
                        Length of k-mer to check after the primer sequence.
                        Defaults to 30.
  --kmer_threshold KMER_THRESHOLD
                        Number of single nucleotide variants deemed acceptable
                        in kmer.
  --primer_bases PRIMER_BASES
                        Number of bases from primer region to use in gapped
                        alignment.Helps with variant calling near the edges of
                        a block. Defaults to 5.
  --proportionthresh N  Keep variants which appear in this proportion of the
                        read pairs. For a given target region, and bin
                        otherwise. Defaults to 0.05.
  --absthresh N         Only keep variants which appear in at least this many
                        read pairs. Defaults to 2.
  --qualthresh N        Minimum base quality score (phred).
  --overlap OVERLAP     Minimum proportion of block which must be overlapped
                        by a read. Defaults to 0.9.
  --max_variants N      Ignore reads with greater than this many variants
                        observed. Defaults to 20.
  --reference FILE      Reference sequences in Fasta format.
  --id_info FILE        File containing rs ID information in VCF format.
  --out FILE            Name of output file containing called variants.
  --log FILE            Logs progress in specified file, defaults to stdout.
  --coverdir COVERDIR   Directory to write coverage files, defaults to current
                        working directory.
  --fast                Use gapped alignment less often, leading to faster run
                        time.
  --snvthresh N         Distance between two single nucleotide variants before
                        going to a gapped alignment.
  --genotype            Compute genotypes for SNVs. Defaults to False.
  --ploidy {1,2}        Ploidy for genotyping 1 = haploid, 2 = diploid.
                        Defaults to 2.
  --error ERROR         Expected base read error rate. Defaults to 0.002.
```
