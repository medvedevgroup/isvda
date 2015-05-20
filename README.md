# ISVDA
Iterative Small Variant Discovery Algorithm
# Prerequisites
- Python 2.7
- g++
- SAMtools
- BWA/Bowtie2 (ISVDA supports BWA and Bowtie2, according to our experiment, BWA runs faster, so the default support is BWA, you can use -e to indicate using Bowtie2)
- Freebayes (To use multi-thread mode of ISVDA, parallel version of Freebayes must be downloaded from: https://github.com/ekg/freebayes/blob/master/scripts/freebayes-parallel, and using --parallel_freebayes to indicate the directory)

# Installation
1. Download repositories
2. cd isvda
3. make

# Usage:
**Quick Usage :**
(This requires samtools, bwa and freebayes in PATH, i.e. you can directly run samtools, bwa and freebayes.)

```
python isvda.py -r whole_genome.fa -p 1.fastq 2.fastq -i 6 -w workspace_directory -t 8 -g --parallel_freebayes=parallel_freebayes_directory
```

**Detail Usage:**

```
python isvda.py -h
```
