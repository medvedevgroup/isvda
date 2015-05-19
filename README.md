# ISVDA
Iterative Small Variant Discovery Algorithm
# Prerequisites
- Python 2.7
- g++
- SAMtools
- BWA/Bowtie2
- Freebayes

# Installation
1. Download repositories
2. cd isvda
3. make

# Usage:
**Quick Usage :**
(This requires samtools, bwa and freebayes in PATH, i.e. you can directly run samtools, bwa and freebayes.)

```
python isvda.py -r whole_genome.fa -p 1.fastq 2.fastq -i 6 -w workspace_directory -t 8 -g
```

**Detail Usage:**

```
python isvda.py -h
```
