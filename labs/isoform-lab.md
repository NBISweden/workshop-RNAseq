---
layout: default
title:  'Isoform Reference'
---

# Reference guided assembly using Cufflinks or StringTie

Assembly of genes and isoforms using Cufflinks and StringTie is a two-step process. 
First you map all the reads from your experiment to the reference sequence. For information on how to map reads to a reference, go to the [mapping reads tutorial](mapping_reads). Then you run another step where you use the mapped reads to assemble potential transcripts and identify the genomic locations of introns and exons.  If you have multipe samples, you can merge the results into one file using Cuffmerge (part of the Cufflinks suite).


## Data available for exercise

In order to make the steps run quickly during the lab, we have extracted only those RNA-seq reads that mapped to one gene (*RAB11FIP5*)
that we will examine in more detail in a later exercise. 

There are files that have been pre-mapped using TopHat and those can be found in

	/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/refBasedAssembly/RAB11FIP5/BAMfiles

on UPPMAX.

There are files that have been pre-mapped and pre-analyzed using TopHat and Cufflinks and those can be found in

	/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/refBasedAssembly/RAB11FIP5/GTFfiles
	
on UPPMAX.
 

## Reference guided assembly using Cufflinks

Cufflinks can do reference-guided assembly, which means that it tries to discover transcripts based on reads mapped to the genome, without considering previous gene annotation (actually there is an option to use annotation as well but we will ignore that for now).
This functionality works even on our small files.

Try to do a reference-guided assembly. This is done simply by running Cufflinks 
without feeding it a GTF file with the -G flag:

     module load cufflinks/2.2.1
     cufflinks -o /path/to/outputDirectory sample.sorted.bam

Substitute the appropriate names for the BAM file and the output directory. When 
Cufflinks has finished (which should hardly take any time at all), the output 
directory will contain a file called ``transcripts.gtf``. Rename the file to a 
name that reflects what created the GTF file.  You can import the GTF file into
IGV as a track for visualization.

&#10067; *Was Cufflinks able to assemble your alignments into something that makes sense?*
 
 
## Reference guided assembly using StringTie

StringTie can also do reference-guided assembly similar to Cufflinks. In their own paper, they claim they can do better assembly than Cufflinks in much shorter time. We have not done any comparison on this ourselves and can not tell if this is the truth. To read the paper go [here](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3122.html). To view a short example on their webpage for examples when it does better go [here](http://ccb.jhu.edu/software/stringtie/#tab3). This functionality works even on our small files.  

Try to do a reference-guided assembly. This is done simply by running StringTie 
without feeding it a GTF file with the -G flag

	module load StringTie/1.3.3
	stringtie -o /path/to/outputDirectory/sample.gtf sample.sorted.bam

Substitute the appropriate names for the BAM file and the output file.
Remember to choose an output file name that reflects what data and programs were used to create it.
When StringTie has finished (which should hardly take any time at all), have a look at the output file.
You can import it into IGV as a track for visualization.

&#10067; *Was StringTie able to assemble your alignments into something that makes sense?*
 

## Merging multiple assemblies into one assembly using Cuffmerge

Cufflinks can merge multiple assemblies into one using a program called cuffmerge. For more information regarding how to run cuffmerge, go [here](http://cole-trapnell-lab.github.io/cufflinks/cuffmerge/). 
This step is not mandatory, but if you want you can try it out. If you have not already loaded the cufflinks module, remember to do so. 

    module load cufflinks/2.2.1
	cuffmerge assembly_GTF_list.txt
	
where ``assembly_GTF_list.txt`` is a text file with a list (one per line) of GTF files that you’d like to merge together into a single GTF file.


## Quantification with Cufflinks

Cufflinks is a well-known software package for estimating gene and isoform abundances in BAM or SAM files based on an annotation file in GTF format. However, we run into problems here because of the small file size. Cufflinks needs a certain amount of data to be able to do its estimations so it will not work very well on our small alignment files. Therefore we have run it on the full alignment files on each sample and provided the merged results at ``/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/isoform_fpkm_table.txt``
(isoform FPKM) and ``/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/RNAseqfpkm_table.txt``.
If you are interested, you can check in the isoform FPKM file how much each isoform 
is deemed to have been expressed by Cufflinks, but we will not go into that today. 
(These results will be revisited tomorrow, in the differential expression lab exercise.)

For reference, the commands we used were of the form

     # Only for reference, does not need to be executed during the exercise
     cufflinks -p 8 -G /proj/uppstore2017171/courses/RNAseqWorkshop/downloads/isoform/otherData/Homo_sapiens.GRCh38.77.fixed.gtf \
	   -o cufflinks_out_137_1 accepted_hits_137_1.bam

The ``-G`` option points to an annotation file in GTF format for which to calculate
FPKM values. The input here is a BAM file (i.e. a binary version of a SAM file) containing RNA-seq read alignments to human genome assembly GRCh38.  

If you were to use your own created GTF file as reference, just replace ``Homo_sapiens.GRCh38.77.fixed.gtf`` with ``merged.gtf`` and run the program:

     # Only for reference, does not need to be executed during the exercise
     cufflinks -p 8 -G /path/to/your/merged/assembly/merged.gtf -o cufflinks_out sample.bam

Other options for doing abundance estimation are [RSEM](http://deweylab.biostat.wisc.edu/rsem/) or the flexible [RPKMforgenes.py script](http://sandberg.cmb.ki.se/media/data/rnaseq/instructions-rpkmforgenes.html).








