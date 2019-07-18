# Run-through

This run-through will show go through the steps of the genomic pipeline (from CBCLs to VCF) using non-DRAGEN tools.

### Computer Specs

| | |
|-|-|
| CPU | Intel Xeon E5-2678 v3 @ 2.50 GHz |
| Cores | 48 |
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
```basespace_NA12878``` is the directory that contains the run data.

```--no-lane-splitting``` to prevent FASTQ files from being split by lane.

```
bcl2fastq -R basespace_NA12878 -o /mnt/genomics/fastq_files_NA12878 --no-lane-splitting
```

__Time Output:__
```42171.47s user 300.86s system 4458% cpu 15:52.54 total```

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
