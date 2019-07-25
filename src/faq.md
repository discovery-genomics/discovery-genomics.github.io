## How to detect pair-ended vs non-pair-ended FASTQ
If pair-ended, there should be 2 files.
	- R1.fastq
	- R2.fastq
The first line represents an id of the read. Information to detect pair-end is at the end of this id. Sequence identifiers format vary by organization. Here is an example of a pair of reads with Illumina sequence identifiers v1.4 and above.

**FT-SA14287_S90_L004_R1_001.fastq**
Notice the 1 near the end of the first line. **1**:N:0:GCCAATAT
```
@A00284:43:HC5TNDSXX:4:1101:1018:1000 1:N:0:GCCAATAT
NTGAAAAGGCAGAATGTATAGTGTTGGTCCTTTGATGGAATTTCTGAGAAAAAACAAACTATTCCAGATGATGACTAAATCACATAGTTTTAAATCTTCTGTAGATTTTTAATGGTTTTTTTTCAAAAACGCATATTGTTCAATATAATAA
+
!FF,FFFFFFFFFFFFFFFF::FF,F:FFFFFFFFFFFF:FFF,FFFFFFFF:FFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFF:FFFFFFFFFFFFFFF:F,FFF:,FFFFF:FFFFFFFFFFFFFFFF
```

  **FT-SA14287_S90_L004_R2_001.fastq**
  Notice the 2 near the end of the first line. **2**:N:0:GCCAATAT
  ```
  @A00284:43:HC5TNDSXX:4:1101:1018:1000 2:N:0:GCCAATAT
TATAAATAAACTTTCTGAAATTTTGTGGTTCAAAAGGATGAATGGTGGGGATTATCATAATCTTAAATTTCTACAAAACTCTAGCAAAGATACTAACTGGTGTGATTTTTTTCTTTTTTGTATCATACTGCACATTTGTGAATCCACCAAG
+
FFFF,FF,FF,FFFF,F:FFF:,FFF,FFFFFFFF::FFFFFFF:FFFFFFFFFFFFFFFFF:FFF::FF:FFFFFFF,F::,::FFFFF:FF,FF,FFFFFFFFFFFF,:FF:,F:FFFFF,:FFF:FF,FFFFF:,FFFFFF,F:FF,:
  ```

Note that the identifier in the first lines must be the same for both files to be pair-ended.

Sometimes pair-ended FASTQ can be in a single file. This is means that the sequence is paired, but interleaved. If there are two reads with the same id, the FASTQ can be paired. Check if there are duplicate ids to verify if the FASTQ is interleaved or not.

Here is an example of a pair of reads with Illumina sequence identifiers using an older version of the format. v1.4 and below.
Note the **/1** and **/2** and the end.

**R1.fastq**
```
@HWUSI-EAS100R:6:73:941:1973#0/1
```
**R2.fastq**
```
@HWUSI-EAS100R:6:73:941:1973#0/2
```
References
- [https://www.biostars.org/p/95310/](https://www.biostars.org/p/95310/)
- [https://en.wikipedia.org/wiki/FASTQ_format](https://en.wikipedia.org/wiki/FASTQ_format)

## Should I use GRCh38 or GRCh38p13?
Use GRCh38p13, the p13 means it is the 13th patched version of GRCh38. Therefore it is a more up to date/accurate version.

## What are read groups?
A [read group](https://software.broadinstitute.org/gatk/documentation/article.php?id=6472) is defined as a set of reads that a sequencer generates in a single run. A read group is effectively treated like a separate run of the instrument in downstream analysis, such as in BQSR, where same read groups are expected to share the same error model.

Here some example:

1. 1 sample that is run on a single lane of a flow cell. All the reads from that lane are part of the same read group. 1 SAM file, with 1 @RG in the header that specifies the read group. All the read records in this file are part of that read group.

2. 4 samples. 1 dad, 1 mom, 1 son, 1 daughter. Each is ran on the same flowcell, but different lanes. Dad on lane1, mom on lane2, etc... This means I will have 4 different SAM files with the corresponding read group in the header.

3. 1 sample, however it is prepped with two different DNA libraries. One with 400 bp inserts and one with 200 bp inserts. Each library is run on two lanes. In this situation, there would be 4 SAM files with different read groups. 400 bp lane 1, 400 bp lane 2, 200 bp lane 3, 200 bp lane 4.

__GATK Required Read Groups__
- __ID: Read Group Identifier__ Identifier that assigns reads to a read group. Each read group's ID must be unique. @RG is in the header, which defines the read group. Each read record will have a RG:Z tag that references the read group ID it belongs to. Illumina read group IDs are in the following format. flowcell + lane name and number.
- __PU: Platform Unit__ Holds 3 types of information, {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. FLOWCELL_BARCODE is the unique id for a flow cell. LANE means which lane within that flowcell. SAMPLE_BARCODE is the sample/library identifier. PU is not required by GATK, but takes precedence over ID for base recalibration.
- __SM: Sample__ Name of the sample in this read group. All read groups with the same SM value are considered as sequencing data from the same sample. If sequencing a pool of samples, use the name of the pool. If multiplexing, keep the SM tag as each sample name. The difference between pooling and multiplexing is that pooling is for samples that are not barcoded.
- __PL: Platform Sequencing Technology__ What sequencing technology was used. Valid values are: ILLUMINA, SOLID, LS454, HELICOS, PACBIO.
- __LB: DNA Prep Library Identifier__ What DNA library did you use.
