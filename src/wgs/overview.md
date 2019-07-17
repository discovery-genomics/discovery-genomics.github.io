# Whole Genome Sequencing

A brief overview of the WGS (Whole Genome Sequence) production and analysis process is
1. Load target DNA into Illumina NextSeq
2. Illumina NextSeq produces raw data in __CBCL__ format
3. Raw __CBCL__ data is processed and turned into industry standard __FASTQ__ format
  - Quality control is performed after this step
4. __FASTQ__ file is preprocessed to remove abnormalities to produce a __Trimmed FASTQ__
5. __Trimmed FASTQ__ is mapped and aligned to reference genome (often a [Genome Reference Consortium Reference](https://www.ncbi.nlm.nih.gov/grc/data)) to produce a __BAM__ file.
6. The __BAM__ is then passed through 2 processes to sort and mark duplicate resulting in a __Processed BAM__.
  - Quality control is performed after this step
7. __Processed BAM__ is then checked for variants using one of several variant calling tools to produce a __VCF__.
8. (If not already recalibrated) The __VCF__ is recalibrated to produce more accurate quality scores.
9. The VCF is ready to be annotated.
