---
layout: default
title:  'Index'
---

## RNA-Seq Lab

This page contains links to different tutorials that are used in the RNA-seq course. Some of the tutorials are well documented and should be easy to follow. We also supply more beta versions of labs that requires more input from the user and may contain errors. 

### 1. Introduction

In the links below, there are information about tools and data that we will used during the other labs. Please make sure you know the data and how to use **R** and **IGV** before you proceed with the other labs. 

*	[Introduction to the RNA-Seq data provided](intro)  
*	[Short introduction  to R](R_intro)  
*	[Short introduction to IGV](IGV) 

### 2. Mapping reads 

This section contains information on how to map reads to a reference genome. In the tutorial you will learn how to use both **STAR** and **HISAT2**.
 
*	[Mapping reads to a reference and converting the results to BAM format](mapping_reads) 

### 3. Transcript assembly

This section contains information regarding how to assemble short reads into transcripts.

*	[Reference guided assembly](isoform-lab)  
*	[*De novo* assembly](isoform-denovo)

### 4. Visualise mapped reads and assembled transcripts on reference

When reads have been mapped to a reference and/or assembled to transcripts, it is always a good idea to check on a reference what the results look like.
 
*	[Isoform visualisation](isoform-visualisation)  

### 5. Quality control laboratory

Before doing any other analysis on mapped RNA-seq reads it is always important to do quality control of your mapped reads and that you do not have any obvious errors in your RNA-seq data. 

*	[RNA-Seq Quality Control](QC_lab)   

### 6. Small RNA analysis

When working with small RNA RNA-seq reads, in this case; miRNA, there are some parts of the analysis that are different. This will be covered in this lab.  

*	[small RNA analysis](smallRNA-lab)

### 7. Differential expression analysis

There are many software packages for differential expression analysis of RNA-seq data.

Several tools, such as DESeq and edgeR, start from read counts per gene and use the discrete nature of the data to do statistical tests appropriate for this kind of data. It can be argued that that such counts will never give quite the correct results because the presence of alernative isoforms confound the read counting. Cuffdiff therefore combines isoform-level quantification and differential expression testing into one framework and claim that they achieve better results because they are able to take into account the uncertainty of isoform quantification. 

*	[Differential expression analysis using DEseq2](DEseq2)
*	[Differential expression analysis using Sleuth](kallisto)
*	[Differential expression analysis using multivariate analysis in SIMCA](Simca_tutorial)

### 8. Gene set analysis

We will perform gene-set analysis on the output from the tutorial on Differential expression analysis of RNA-seq data using DESeq.

*	[Gene set analysis](GSA_tutorial)  

### 9. Beta labs 

There are some labs that are more close to the cutting edge of analysis and therefore are not as well tested as the ones above. These are tools that have high potential and will most likely, if they hold, will be moved to the mature labs.
 
*	[Single cell RNA PCA and clustering](PCA_clustering_single_cell)    

 
### 10. UPPMAX
 
 One example of a sbatch script
 
 *  [sbatch scripts](sbatchScript)   
  
 
### 11. Caveat

We will try to keep these tutorials up to date. If you find any errors, or things that you think should be updated, please contact Johan (johan.reimegard@scilifelab.se).
  		
