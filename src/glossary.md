## File Formats
**BCL:** Binary Call Format.

**CBCL:** Concatenated Binary Call Format.

**FASTA:** File format used to store reference genome.

[**FAI:**](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference) FASTA Index File. Used for fast random access. Needed in order to use the FASTA file as reference in GATK HaplotypeCaller.

[**DICT:**](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference) Dictionary file. Used for fast random access. Needed in order to use the FASTA file as reference in GATK HaplotypeCaller.

**FASTQ:** Stores genetic sequence fragments with quality scores.

**SAM:** Sequence Alignment Map. Has information about each read and the reference sequence it maps to. Human readable.

**BAM:** If aligned, then it is a compressed SAM. If not aligned, then it is more like a compressed FASTQ

[**BAI:**](https://www.biostars.org/p/15847/) Companion index file for a BAM. Acts like a table of contents for a BAM, which allows programs to jump to an index in the respective BAM file. Needed for variant calling with GATK HaplotypeCaller.

**uBAM**: Unaligned BAM. Compressed FASTQ

**VCF:** Variant Call Format. Contains the reads that differ from the reference genome.

**VCF.GZ.TBI:** Index file of vcf.gz

[**GVCF:**](https://gatkforums.broadinstitute.org/gatk/discussion/4017/what-is-a-gvcf-and-how-is-it-different-from-a-regular-vcf) VCF with extra information.

**FNA:** FastA format file containing Nucleotide sequence (DNA)

**GBFF:** Genbank Genome file containing genome sequence and annotation

**GFF:** general feature format containing genomic regions, the "genes, transcripts, etc"

**FAA:** FastA format file containing Amino-acid sequence (Protein, peptide)

**GPFF:** Genbank Protein file containing protein sequence and annotation

**\*.gz:** Compressed file. Some tools can work off of compressed files and output compressed files.
## Data Conversion Tools
[**bcl2fastq**](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf)
Illumina software that demultiplexes raw sequence data (BCL) to FASTQ.

[**bwa**](http://bio-bwa.sourceforge.net/bwa.shtml)
Burrow-Wheeler Aligner.
A program for mapping/aligning sequences with a reference genome.

[**samtools**](http://www.htslib.org/doc/samtools.html)
Program that provides various utilities to manipulate SAM and BAM files. Some features include converting between SAM and BAM file formats, sorting, indexing, filtering.

[**picard**](https://broadinstitute.github.io/picard/)
Like samtools, picard provides a set of tools to manipulate SAM, BAM, CRAM, and VCF files. Uses Java 8.
[more info about picard and metrics](https://www.broadinstitute.org/files/shared/mpg/plathumgen/plathumgen_fennell.pdf)

[**gatk**](https://software.broadinstitute.org/gatk/)
Toolkit that provides many genomic analysis tools. Mainly focused on variant calling and genotyping.
- [**HaplotypeCaller**](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php)
Makes variant calls on a BAM file. The VCF that is outputted must be filtered by variant recalibration or hard-filtered before use in down-stream analysis. Can also output a GVCF, in which case it must be ran through GenotypeGVCFs and then filtering.
Requires the FASTA file to be indexed with a [FAI and DICT file.](https://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference)
Requires the BAM file to be indexed with a BAI file (samtools index).

**mrsfast**


preprocessing
trimming fastq

sort
index

deduplicate/mark duplicate

## Run Data Files and Directories

```
basespace_NA12878
├── Config
├── Data
├── InterOp
├── NovaSeq 6000 SP: TruSeq PCR-Free 450 (2 replicates of NA12878)_166552386.json
├── Recipe
├── RTA3.cfg
├── RTAComplete.txt
├── RunInfo.xml
├── RunParameters.xml
├── SampleSheet.csv
└── SequenceComplete.txt
```

cbcl files are in ```Data/Intensities/BaseCalls```

```SampleSheet.csv``` contains information about samples, library chemistry, adapter trimming sequences, and other information about the run. ```bcl2fastq``` refers to this file to demultiplex BCL to FASTQ. For most runs, a sample sheet is optional. If no sample sheet is provided, all reads are written to a FASTQ file named similar to ```Undertermined_S0```, which is the default sample for reads that are not assignable a sample (ex: poor quality, wrong indexes in sample sheet, sequencing error).
