---
layout: default
title:  'Isoform Visualisation'
---

# Visualizing isoforms in RNA-seq data


One of the advantages of RNA-seq data compared to microarrays is that you get 
isoform-level information 'for free' in addition to gene-level information. 
In this exercise, we are going to look at some ways to identify and visualize isoforms.

 As you have hopefully read on the introduction page, we will use RNA-seq and quantitative 
mass-spectrometry (MS) data sets from the A431 cell line. These data were measured by a 
group at SciLifeLab Stockholm in a 'proteogenomics' study where the aim was to discover 
new genes or gene variants by deep proteomic profiling using an MS method, and mapping 
the obtained peptides back to the genome. 
The RNA-seq data was obtained to see if there was RNA-level support for the predicted novel 
peptides and transcript variants. We will look at one locus that was flagged by the research 
group as being interesting, and see what the RNA-seq data look like for that gene.


## Strategies for using the RNA-seq data

There are different ways to find out how the RNA-seq data shows the *RAB11FIP5* gene to 
be expressed. Roughly speaking, we can talk about three different strategies:

*	Mapping the reads to a reference genome and quantifying the gene and isoform FPKMs  
*	Mapping the reads to a reference genome and assembling transcripts based on the mapping (reference guided assembly)  
*	Assembling the reads *de novo*, ignoring the reference genome for now  


In order to make these steps computationally feasible during the lab, we have extracted 
only those sequences that mapped to the *RAB11FIP5* gene in each sample. These "sub-FASTQ" 
files can be found in ``/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/RAB11FIP5_fastqFiles``.
 

To do the reference guided assembly yourself go to [Reference guided assembly using Cufflinks or Stringtie](isoform-lab). 

This link also contains information on how to quantify already annotated genes and isoforms.

To do the *de novo* assembly yourself go to [Isoform detection using RNA-seq *de novo* Assembly](isoform-denovo).



## Visualise isoform data

For the gene *RAB11FIP5* you will now hopefully have generated your own data that you can look at. 
If everything worked you will now have:

 * One or two BAM files with short reads mapped to the genome using HISAT2 and/or STAR 

 * One GTF file containing different isoforms of *RAB11FIP5* assembled with Cufflinks and/or Stringtie based on the short reads that mapped to the genome
 
 * One BAM file with **Trinity** assembled transcripts mapped to the genome


### Importing reference based isoform info for the *RAB11FIP5* gene

Since it takes time to generate all data, we have already created other files that you can also download and view in your browser. This includes result files for the subset of reads that map to the *RAB11FIP5* gene. These mappings have been used for reference based assembly of isoforms. 

You can find all BAM files and GTF files for all samples here  `/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/refBasedAssembly/RAB11FIP`. You can then view them in IGV using "Import from URL".

### Importing de novo assembled transcripts mapped to the *RAB11FIP5* gene

We have also created result files from *de novo* transcriptome assembly using the the subset of reads that map to the *RAB11FIP5* gene. The assembled transcripts were
then mapped back to the genome. 

You can download all BAM files and GTF files for all samples from uppmax from here `/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/deNovo/BAMfiles`


## Importing reference based isoform info to the genome

In addition, we provide result files for *all* reads that were mapped to the genome. These mappings were used for 
reference based assembly of isoforms across the entire genome. There is a GTF file with the final merged isoform  
information from all 12 samples.  

You can download all BAM files and GTF files for all samples from uppmax from here `/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/refBasedAssembly/Genome`. 

## Importing the peptide track for the *RAB11FIP5* gene and the genome                                                           

As mentioned above, we will look at some identified peptides from a mass-spectrometry 
experiment, and compare those with RNA-seq data from the same cell line.

You can download the BED file containing all peptides mapped to the genome from uppmax `/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData` .


From the name of the BED file (human_A431_global-TDA-FDR1pc_green-known_red-novel.bed), IGV will automatically know to color the track according to peptide status
(green for annotated peptides, red for novel peptides).


# Importing the Pac bio reads mapped to the genome                                                         

You can use a web browser to access the BAM file containing all PacBio reads mapped to the genome  from uppmax `/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/).

How do the PacBio reads from full length transcripts compare to the transcripts assembled from short reads?






















