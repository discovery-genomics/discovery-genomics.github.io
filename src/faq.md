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
