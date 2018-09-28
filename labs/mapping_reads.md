---
layout: default
title:  'Mapping reads'
---


# Mapping reads to a reference and converting the results to BAM format


## Data available for exercise

You can carry out this exercise using the RNA-seq data we provide or your own data, if you have any.

In order to make the steps run quickly during the lab, we have extracted only those RNA-seq reads that mapped to one gene (*RAB11FIP5*)
that we will examine in more detail in later exercises. 

All FASTQ files that you will need for this exercise can be found in
 
	/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/referenceBased/data

on UPPMAX.

If you want to map more files for practice, you can continue with files for additional samples, found in
	
	/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/RAB11FIP5_fastqFiles

on UPPMAX.

A pre-built human genome index for HISAT2 is found here
 
	/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/reference/hg19_hisat2

on UPPMAX.
Note that *hg19* indicates the version of the human genome assembly that was indexed.

A pre-built human genome index for STAR is found here
 
	/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/reference/hg19_Gencode14.overhang75

on UPPMAX.
 

## Mapping short reads to a reference using HISAT2

Here, you will map the reads to the hg19 reference genome using the RNA-seq aligner **HISAT2**. Note that if you are using your own non-human data, you need to use a reference genome for the corresponding species.

There are many features that can be tweaked using HISAT2. For more information on all flags that can be used go [here](https://ccb.jhu.edu/software/hisat2/manual.shtml).

Read below for the flags we use for this exercise. Remember to change filenames accordingly so that you can run your program properly and know which files you have used.

To load the HISAT2 module on UPPMAX, execute:
     
     module load bioinfo-tools      # This is to get access to all bioinformatics tools available on UPPMAX
     module load HISAT2/2.1.0       # The specific mapping program

Now you can map the reads from one of the samples (or several; it's up to you which ones) using a command such as the one below.

	mkdir outDir
    
    hisat2 -p N --dta-cufflinks -x path/to/index/fileName \
	  -1 path/to/reads/sample_1.fastq -2 path/to/reads/sample_2.fastq \
	  -S outDir/hisat2.sam --summary-file outDir/hisat2_summary.txt
    
The flags used are:

*  ``-p N`` specifices the number of threads that will be used by the program
*  ``--dta-cufflinks`` will generate output that is optimal for downstream analysis with Cufflinks
* ``-x /path/to/index/fileName`` specifices the path to the pre-built genome index. Note that the index consists of multiple files ending in ``.ht2``, and only the shared part of the filename should be indicated (e.g. ``genome`` if the files are called ``genome.1.ht2``, ``genome.2.ht2``, etc).
*  `` -1 /path/to/reads/sample_1.fastq `` is where you should list the first-read FASTQ files that you wish to map 
*  `` -2 /path/to/reads/sample_2.fastq `` is where you should list the second-read FASTQ files that you wish to map
*  ``-S outDir/hisat2.sam`` is the name of the result file that will be created
*  ``--summary-file outDir/hisat2_summary.txt`` is the name of a file for summary information about the alignments

This should run fairly quickly and create the files you specified with ``-S`` and ``--summary-file``.

If everything worked, HISAT2 should report some statistics about how many reads were mapped, on your terminal and in the summary file.

&#10067; *Try to answer the following:*

* How many RNA-seq read pairs were provided as input to HISAT2?
* How many of those read pairs were mapped by HISAT2?
* How many reads were uniquely mapped, i.e. mapped to one genomic location?
* In general, do the alignments seem to be good? I.e. do they cover the entire reads and contain few mismatches?

To answer these questions, you should look at the input to and output from HISAT2. You may also need to consult the [HISAT2 manual](https://ccb.jhu.edu/software/hisat2/manual.shtml), [information about the FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) and the [SAM format specification](https://github.com/samtools/hts-specs).

## Mapping short reads to a reference using STAR

Here we will map the reads to the hg19 reference genome using a popular RNA-seq 
aligner, **STAR**. There are many many features that can be tweaked using STAR. For more information concerning different features that can be used see the [manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf).

Read below for the flags we use for this exercise. Remember to change filenames accordingly 
so that you can run your program properly and know which files you have used.

To load the STAR module on UPPMAX, execute

     module load bioinfo-tools
     module load star
     
  
Now you can map the reads from one of the samples (or several; it's up to you 
which ones) using a command such as the one below.
  
	mkdir outDir
    
	STAR  --runThreadN N --outSAMstrandField intronMotif --genomeDir /path/to/index \
	  --readFilesIn /path/to/reads/sample_1.fastq /path/to/reads/sample_2.fastq \
	  --outFileNamePrefix outDir/
	
flags used are 

*  ``--runThreadN N`` specifies the number of threads that will be used by the program.
*  ``--outSAMstrandField intronMotif`` adds information  (to the SAM output file) required for downstream analysis with Cufflinks
*  ``--genomeDir /path/to/index`` specifies the directory containing the pre-built genome index
*  ``--readFilesIn /path/to/reads/sample_1.fastq /path/to/reads/sample_2.fastq`` is where you should list the FASTQ files that you wish to map
*  ``--outFileNamePrefix outDir`` specifies the output directory

This should run fairly quickly and put a file called ``Aligned.out.sam`` in 
the directory that you specified with ``--outFileNamePrefix``. 

&#10067; *Look at the output files from STAR, and try to answer the same questions as for HISAT2 above.*

## Converting SAM files to BAM files

If you were able to run HISAT2 and STAR sucessfully, this should have produced files with mapped reads in SAM format. These files need to be converted to *sorted* and *indexed* BAM files for efficient downstream analysis.

You should try to give the BAM files representable names, in order to make it easier to manage your files. A good naming scheme for BAM files is to use names that indicate what you mapped and how. As an example, if you mapped sample 12 using HISAT2 you could create a file named ``sample12_RAB11FIP5.hg19.HISAT2.bam``. 

The most commonly used tool for converting from SAM to BAM is [Samtools](http://www.htslib.org/doc/samtools.html) (follow the link for more information about Samtools).

To load the Samtools module on UPPMAX, execute:

    module load bioinfo-tools
    module load samtools

The Samtools command to convert from SAM to a sorted BAM file is:

	samtools sort -o output.bam input.sam

Remember to use an appropriate filename instead of ``output.bam``!

Next, we need to index the BAM file.

	samtools index properName.bam

The indexing step should create an index file with the suffix ``.bai``.

You can also get a report on your mapped reads using the samtools command *flagstat*:

	samtools flagstat properName.sorted.bam > properName.sorted.flagstat

Since the BAM file contains all the information from the original SAM file, remember to remove the SAM file once you are finished, in order to free up disk space.

The sorted, indexed bam file can be viewed in the Integrative Genomics Viewer (IGV). Instructions are [here](IGV).

&#10067; *Try to answer the following:*

* Does the output from ``samtools flagstat`` confirm any of your answers to the questions in the HISAT2 and STAR sections above?
* Load the the BAM files with HISAT2 and STAR results into IGV. Go to the RAB11FIP5 locus. Have HISAT2 and STAR mapped the reads in a similar way?
* Detailed examination of the read alignments in IGV should indicate if the RNA-seq data is strand-specific. Is it?

## When you are done...

...compare your answers to the solutions [here](mapping_reads_answers).
