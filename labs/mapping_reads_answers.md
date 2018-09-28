---
layout: default
title:  'Mapping reads'
---

# Read mapping: Answers to exercise questions

## Mapping short reads to a reference using HISAT2

How many RNA-seq read pairs were provided as input to HISAT2?

	768 pairs (768*2 = 1536 reads)
	From HISAT2 summary output file, or input FASTQ file.

How many of those read pairs were mapped by HISAT2?

	693, sum of:
	  682 concordantly exactly 1 time
      0 concordantly >1 times (i.e. multimapping)
      11 discordantly
	(from HISAT2 summary output)

How many reads were uniquely mapped, i.e. mapped to one genomic location?

	1440
	
	Can compute this by running:
	  grep -cw 'NH:i:1' hisat2.sam
	
    Or if you have made a BAM file:
      samtools view hisat2.bam | grep -cw 'NH:i:1'
	
    Or by looking at the summary output file, summing up:
      2 * 682 concordantly and uniquely mapped pairs
      2 * 11 discordantly and uniquely mapped pairs
      54 uniquely mapped reads lacking a mapping for the paired read

In general, do the alignments seem to be good? I.e. do they cover the entire reads and contain few mismatches?

    Yes.
	
	Look at mapping quality, CIGAR string and NM tag (edit distance) in the SAM file:
     cut -f5-6,12- hisat2.sam | more
	
	Or for the BAM file:
     samtools view hisat2.bam | cut -f5-6,12- | more

## Mapping short reads to a reference using STAR

How many RNA-seq read pairs were provided as input to STAR?

	768 pairs (768*2 = 1536 reads)
	(From Log.final.out, or input FASTQ file)

How many of those read pairs were mapped by STAR?

	738, sum of
	  733 uniquely mapped
	  5 multimapped
	(From Log.final.out)

How many reads were uniquely mapped, i.e. mapped to one genomic location?

	1466
	
	This can be obtained from Log.final.out (2 * 733 uniqely mapped reads).

	Or by running:
	  grep -cw 'NH:i:1' Aligned.out.sam
	
	Or on a BAM file:
	  samtools view star.bam | grep -cw 'NH:i:1'

In general, do the alignments seem to be good? I.e. do they cover the entire reads and contain few mismatches?

	Yes.
	See Log.final.out, in particular "Average mapped length" and "Mismatch rate per base".
	Can also use samtools (see the answer to the corresponding question for HISAT above).

## Converting SAM files to BAM files

Does the output from samtools flagstat confirm any of your answers to the questions in the HISAT2 and STAR sections above?

	Kind of. Note that several of the numbers refer to the number of mappings (not reads) and there can be several mappings per read.

Load the the BAM files with HISAT2 and STAR results into IGV. Go to the RAB11FIP5 locus. Have HISAT2 and STAR mapped the reads in a similar way?

	Yes.

Detailed examination of the read alignments in IGV should indicate if the RNA-seq data is strand-specific. Is it?

	No.
	Right-click on the track name, choose "Color alignments by" -> "First of pair strand".
	If the data is strand-specific, one color should dominate for every gene (depending on the strand of the gene).
