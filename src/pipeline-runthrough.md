# Run-through

This run-through will show go through the steps of the genomic pipeline (from CBCLs to VCF) using non-DRAGEN tools.

### Computer Specs

| | |
|-|-|
| CPU | Intel Xeon E5-2678 v3 @ 2.50 GHz |
| Cores | 24 |
| RAM | 126 GB |
| Memory | 37 TB |

### Tools Used
- bcl2fastq

### Initial Files

- [__Reference Genome GRCh38.p13__](https://www.ncbi.nlm.nih.gov/assembly/?term=GRCh38)
- __Run Data__
  - Access the [BaseSpace public data sets](https://basespace.illumina.com/datacentral)
  - Search _NovaSeq 6000 SP: TruSeq PCR-Free 450 (2 replicates of NA12878)_
  - Import the run to your account
  - Download via [BaseSpace CLI](https://developer.basespace.illumina.com/docs/content/documentation/cli/cli-overview)
  - The run data should be about 94 GB

### Convert the run data (CBCLs) to FASTQ
```
bcl2fastq -R basespace_NA12878 -o /mnt/genomics/fastq_files_NA12878 --no-lane-splitting
```
basespace_NA12878 is the directory that contains the run data.

```
--no-lane-splitting
```
to prevent FASTQ files from being split by lane.




__Time Output:__
```
42171.47s user 300.86s system 4458% cpu 15:52.54 total
```

__fastq_files_NA12878 Contents__
```
fastq_files_NA12878
├── NA12878-PCRF450-1_S1_R1_001.fastq.gz
├── NA12878-PCRF450-1_S1_R2_001.fastq.gz
├── NA12878-PCRF450-2_S2_R1_001.fastq.gz
├── NA12878-PCRF450-2_S2_R2_001.fastq.gz
├── Reports
├── Stats
├── Undetermined_S0_R1_001.fastq.gz
└── Undetermined_S0_R2_001.fastq.gz
```

Decompress with gzip to get the FASTQ files.
Four decompressed FASTQ files add up to 656 GB.
Meaning 94 GB of CBCLs converted to 656 GB of FASTQs.

### Map/Align FASTQ to BAM

1. Create bwa index files for the reference FASTA file.
```
bwa index -a bwtsw GCF_000001405.39_GRCh38.p13_genomic.fna
```
  __Time Output__
  ```
  4208.44s user 26.93s system 99% cpu 1:10:47.24 total
  ```

2. Align sample 1 and sample 2 as pair-ended reads, which will create a SAM file (which is too large to store). Pipe the output to samtools to end up with a BAM.

Ran both samples in parallel, half the threads for sample 1 and half the threads for sample 2.

__Sample 1__
```
bwa mem -t 24 GCF_000001405.39_GRCh38.p13_genomic.fna fastq_files_NA12878/NA12878-PCRF450-1_S1_R1_001.fastq fastq_files_NA12878/NA12878-PCRF450-1_S1_R2_001.fastq | samtools view -@ 24 -o NA12878_S1.bam -
```
  __Time Output__
  ```
  bwa mem -t 24 GCF_000001405.39_GRCh38.p13_genomic.fna    870775.86s user 6008.81s system 2314% cpu 10:31:23.95 total
  samtools view -@ 24 -o NA12878_S1.bam  23229.54s user 975.26s system 63% cpu 10:31:23.77 total      
  ```
__Sample 2__
```
bwa mem -t 24 GCF_000001405.39_GRCh38.p13_genomic.fna fastq_files_NA12878/NA12878-PCRF450-2_S2_R1_001.fastq fastq_files_NA12878/NA12878-PCRF450-2_S2_R2_001.fastq | samtools view -@ 24 -o NA12878_S2.bam -
```
  __Time Output__
  ```
  bwa mem -t 24 GCF_000001405.39_GRCh38.p13_genomic.fna    762165.59s user 4784.47s system 2312% cpu 9:12:51.64 total
  samtools view -@ 24 -o NA12878_S2.bam  19139.70s user 810.07s system 60% cpu 9:12:51.48 total

  ```
Note: running at ~90 GB of memory.
