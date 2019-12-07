---
layout: default
title:  'Mapping reads'
---

# Mapping: Answers

## Mapping using HISAT2

Alignment summary

> `768 reads; of these:`
>  `768 (100.00%) were paired; of these:`
>    `94 (12.24%) aligned concordantly 0 times`
>    `674 (87.76%) aligned concordantly exactly 1 time`
>    `0 (0.00%) aligned concordantly >1 times`
>    `----`
>    `94 pairs aligned concordantly 0 times; of these:`
>      `12 (12.77%) aligned discordantly 1 time`
>    `----`
>    `82 pairs aligned 0 times concordantly or discordantly; of these:`
>      `164 mates make up the pairs; of these:`
>        `104 (63.41%) aligned 0 times`
>        `59 (35.98%) aligned exactly 1 time`
>        `1 (0.61%) aligned >1 times`
> `93.23% overall alignment rate`

How many RNA-seq read pairs were provided as input to HISAT2?

> `768 pairs (768*2 = 1536 reads)`
>
> From HISAT2 summary output file, or input FASTQ file.

How many of those read pairs were mapped by HISAT2?

> `1432 = ((674+12)*2) + (59+1)`

How many reads were uniquely mapped, i.e. mapped to one genomic location?

> `1431`
>
> This can compute this by running:
>
> `grep -cw 'NH:i:1' hisat2.sam`
> 	
> Or if you have made a BAM file:
> 
> `samtools view hisat2.bam | grep -cw 'NH:i:1'`
> 	
> Or by looking at the summary output file, summing up:
> 
> `2 * 674 concordantly and uniquely mapped pairs`
> `2 * 12 discordantly and uniquely mapped pairs`
> `59 uniquely mapped reads lacking a mapping for the paired read`

In general, do the alignments seem to be good? I.e. do they cover the entire reads and contain few mismatches?

> Yes.
	
Look at mapping quality, **CIGAR** string and **NM** tag (edit distance) in the SAM file:
     
> `cut -f5-6,12- hisat2.sam | more`
>	
> Or for the BAM file:
>
> `samtools view hisat2.bam | cut -f5-6,12- | more`

## Mapping using STAR

Alignment summary


>```
>                                 Started job on |       Oct 22 22:45:41
>                             Started mapping on |       Oct 22 22:47:03
>                                    `Finished on |       Oct 22 22:47:06
>       Mapping speed, Million of reads per hour |       0.92
>
>                          Number of input reads |       768
>                      Average input read length |       202
>                                    UNIQUE READS:
>                   Uniquely mapped reads number |       733
>                        Uniquely mapped reads % |       95.44%
>                          Average mapped length |       199.98
>                       Number of splices: Total |       81
>            Number of splices: Annotated (sjdb) |       68
>                       Number of splices: GT/AG |       81
>                       Number of splices: GC/AG |       0
>                       Number of splices: AT/AC |       0
>               Number of splices: Non-canonical |       0
>                      Mismatch rate per base, % |       0.26%
>                         Deletion rate per base |       0.00%
>                        Deletion average length |       0.00
>                        Insertion rate per base |       0.00%
>                       Insertion average length |       1.33
>                             MULTI-MAPPING READS:
>        Number of reads mapped to multiple loci |       5
>             % of reads mapped to multiple loci |       0.65%
>        Number of reads mapped to too many loci |       0
>             % of reads mapped to too many loci |       0.00%
>                                  UNMAPPED READS:
>       % of reads unmapped: too many mismatches |       0.00%
>                 % of reads unmapped: too short |       3.91%
>                     % of reads unmapped: other |       0.00%
>                                  CHIMERIC READS:
>                       Number of chimeric reads |       0
>                            % of chimeric reads |       0.00%
>```

How many RNA-seq read pairs were provided as input to STAR?

> `768 pairs (768*2 = 1536 reads)`
> 
> from **Log.final.out**, or input FASTQ file

How many of those read pairs were mapped by STAR?

> `738, sum of`
> `733 uniquely mapped`
> `5 multimapped`
>
> from **Log.final.out**

How many reads were uniquely mapped, i.e. mapped to one genomic location?

> `1466`
>	
> This can be obtained from **Log.final.out** (2 * 733 uniqely mapped reads).
>
> Or by running:
>	  
> `grep -cw 'NH:i:1' Aligned.out.sam`
>	
> Or on a BAM file:
>
> `samtools view star.bam | grep -cw 'NH:i:1'`

In general, do the alignments seem to be good? I.e. do they cover the entire reads and contain few mismatches?

> Yes.
	
See **Log.final.out**, in particular **Average mapped length** and **Mismatch rate per base**. Can also use samtools (see the answer to the corresponding question for HISAT above).

## Converting SAM to BAM

Does the output from samtools flagstat confirm any of your answers to the questions in the HISAT2 and STAR sections above?

> Kind of. Note that several of the numbers refer to the number of mappings (not reads) and there can be several mappings per read.

Load the the BAM files with HISAT2 and STAR results into IGV. Go to the **RAB11FIP5** locus. Have HISAT2 and STAR mapped the reads in a similar way?

> Yes.

Detailed examination of the read alignments in IGV should indicate if the RNA-seq data is strand-specific. Is it?

> No.
> 
> Right-click on the track name, choose **Color alignments by > First of pair strand**.
If the data is strand-specific, one color should dominate for every gene (depending on the strand of the gene).

