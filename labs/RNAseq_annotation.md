---
layout: default
title:  'Exercise - Functional annotation'
---

# Functional annotation

Functional annotation is the process during which we try to put names to faces - what do the transcript that we have assemble do? Basically all existing approaches accomplish this by means of similarity. If a translation product has strong similarity to a protein that has previously been assigned a function, the function in this newly annotated transcript is probably the same. Of course, this thinking is a bit problematic (where do other functional annotations come from...?) and the method will break down the more distant a newly annotated genome is to existing reference data. A complementary strategy is to scan for more limited similarity - specifically to look for the motifs of functionally characterized protein domains. It doesn't directly tell you what the protein is doing exactly, but it can provide some first indication.

In this exercise we will use an approach that combines the search for full-sequence similarity by means of 'Blast' against large public databases with more targeted characterization of functional elements through the InterproScan pipeline. Interproscan is a meta-search engine that can compare protein queries against numerous databases. The output from Blast and Interproscan can then be used to add some information to our annotation.

There is a pipeline called [trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki) that can now run all the annotation of transcript. Before you had to do each step one by one. We will do a variant of the trinotate pipeline step by step.

## Prepare the input data

For this exercise we will use the transcripts_stringtie.fa output that we create in the assembly part. if you have time, you can also redo the tutorial and annotate the Trinity results (as the assembly was done for the all library and not only for the chr4, it will take a lot more time).

create a new folder for this exercise:
```
cd ~/RNAseq_assembly_annotation/
mkdir RNAseq_annotation
cd RNAseq_annotation
```
You can link the results of Stringtie in your folder

```
mkdir stringtie
cd stringtie/
ln -s ~/RNAseq_assembly_annotation/RNAseq_assembly/stringtie/transcripts_stringtie.fa .
ln -s ~/RNAseq_assembly_annotation/RNAseq_assembly/stringtie/stringtie.gtf .
```

##	Determining longest Open Reading Frames (ORF)

The first step of the annotation of transcript is to determine the longest open reading frame, the longest will be then annotated.
In order to perform this work we use [TransDecoder](https://github.com/TransDecoder/TransDecoder/wiki)
This tool has proved to be worth particularly within the Trinotate pipeline.

TransDecoder identifies likely coding sequences based on the following steps:

1) All open reading frames (ORFs) above a minimum length (default: 100 amino acids) are identified across all transcripts. Incomplete ORFs are accepted (no start or/and no stop)

2) The top 500 longest ORFs are selected and a reading frame-specific
5th-order Markov model is trained based on these coding sequences.

3) All previously identified ORFs are scored as a sum of the log odds ratio across the length of the putative coding sequence. This log-likelihood score is similar to what is computed by the GeneID software.

4) In addition to reporting ORFs meeting the score requirements, any ORF found to exceed a minimum length of 300 amino acids is reported.


```
module load bioinfo-tools
module load TransDecoder
TransDecoder.LongOrfs -t transcripts_stringtie.fa

```
The file that will be used in the rest of the tutorial is in Trinity.fasta.transdecoder_dir/longest_orfs.pep

## Interproscan approach
 Interproscan combines a number of searches for conserved motifs and curated data sets of protein clusters etc. This step may take fairly long time. It is recommended to paralellize it for huge amount of data by doing analysis of chunks of tens or hundreds proteins.

### Perform [InterproScan](https://github.com/ebi-pf-team/interproscan/wiki) analysis
InterproScan can be run through a website or from the command line on a linux server. Here we are interested in the command line approach.
<u>Interproscan allows to look up pathways, families, domains, sites, repeats, structural domains and other sequence features.</u>

Launch Interproscan with the option -h if you want have a look about all the parameters.

- The '-app' option allows defining the database used. Here we will use the PfamA,ProDom,SuperFamily and PIRSF databases.
- Interproscan uses an internal database that related entries in public databases to established GO terms. By running the '-goterms' option, we can add this information to our data set.
- If you enable the InterPro lookup ('-iprlookup'), you can also get the InterPro identifier corresponding to each motif retrieved: for example, the same motif is known as PF01623 in Pfam and as IPR002568 in InterPro.
- The option '-pa' provides mappings from matches to pathway information (MetaCyc,UniPathway,KEGG,Reactome).

```
module load InterProScan

module load perl
interproscan.sh -i transcripts_stringtie.fa.transdecoder_dir/longest_orfs.pep -t p -dp -pa -appl Pfam,ProDom-2006.1,SuperFamily-1.75,PIRSF-3.02 --goterms --iprlookup
```

***Why? What is the error message displaying?***
 <br />

If you did not have a look at the longest_orfs.pep, please have look and find a solution to
make interproscan run.


<br />

<details>
<summary> **Interproscan problem** - Click to expand the solution </summary>

Interproscan is really selective on the fasta input data, there should not be any stop codon * or any character other than ATCG (except in the header of course)

You need to delet the * at in the sequences, you can do :
```
sed -e 's/*//g' transcripts_stringtie.fa.transdecoder_dir/longest_orfs.pep > longest_orfs_stringtie_interpro.fa

```
</details>
 <br />

You can rerun the interproscan with the proper fasta file.

The analysis shoud take 2-3 secs per protein request - depending on how many sequences you have submitted, you can make a fairly deducted guess regarding the running time.
You will obtain 3 result files with the following extension '.gff3', '.tsv' and '.xml'. Explanation of these output are available [>>here<<](https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats).


### load the retrieved functional information in your assembly file:

When you want to use transcript information in an annotation you can either use a fasta file or a gff file. The advantage of the gff file is to be able to merge different functional annotation into it.

Here we are going to create a gff file from the gtf file of stringtie and then add the interproscan annotation to it. We will also reformat the gff so it will be processed properly in the next steps.

```
export PERL5LIB=$PERL5LIB:~/RNAseq_assembly_annotation/GAAS/annotation/
module load BioPerl
~/RNAseq_assembly_annotation/GAAS/annotation/Tools/bin/gxf_to_gff3.pl -g transcripts.gtf -o transcripts.gff3
sed -i 's/transcript/mRNA/g' transcripts.gff3

```

You could write scripts of your own to merge interproscan output into your annotation. But some genome annotation tools (like [Maker](http://gmod.org/wiki/MAKER) in that case) comes with utility scripts that can take InterProscan output and add it to a gff annotation file (you need to load maker).

- ipr\_update\_gff: adds searchable tags to the gene and mRNA features in the GFF3 files.

```
module load maker/3.01.2-beta
ipr_update_gff transcripts.gff3 longest_orfs_stringtie_interpro.fa.tsv > stringtie_with_interpro.gff
```
Where a match is found, the new file will now include features called Dbxref and/or Ontology_term in the gene and transcript feature field (9th column).


## BLAST approach
Blast searches provide an indication about potential homology to known proteins.
A 'full' Blast analysis can run for several days and consume several GB of Ram. Consequently, for a huge amount of data it is recommended to parallelize this step doing analysis of chunks of tens or hundreds proteins. This approach can be used to give a name to the genes and a function to the transcripts.

### Perform Blast searches from the command line on Uppmax:

To run Blast on your data, use the Ncbi Blast+ package against a Drosophila-specific database (included in the folder we have provided for you, under **~/RNAseq_assembly_annotation/assembly_annotation/database/uniprot_dmel/uniprot_dmel.fa**) - of course, any other NCBI database would also work:
```
module load blast/2.7.1+
blastp -db ~/RNAseq_assembly_annotation/assembly_annotation/database/uniprot_dmel/uniprot_dmel.fa -query longest_orfs_stringtie_interpro.fa -outfmt 6 -out blast.out -num_threads 8
```
Against the Drosophila-specific database, the blast search takes about 2 secs per protein request - depending on how many sequences you have submitted, you can make a fairly deducted guess regarding the running time.

### Process the blast outout with Annie
The Blast outputs must be processed to retrieve the information of the closest protein (best e-value) found by Blast. This work will be done using [annie](https://github.com/genomeannotation/annie).

First download annie:
```
cd ~/RNAseq_assembly_annotation/
git clone https://github.com/genomeannotation/Annie.git
```
Then you should load python:
```
module load python/3.5.0
```
Now go back to where you are running the functional annotation analyses and launch annie:
```
~/RNAseq_assembly_annotation/Annie/annie.py -b blast.out -db ~/RNAseq_assembly_annotation/assembly_annotation/database/uniprot_dmel/uniprot_dmel.fa -ipr longest_orfs_stringtie_interpro.fa.tsv -g transcripts.gff3 -o stringtie_annotation.annie
```

Annie writes in a 3-column table format file, providing gene name and mRNA product information. The purpose of annie is relatively simple. It recovers the information in the sequence header of the uniprot fasta file, from the best sequence found by Blast (the lowest e-value).

### load the retrieved information in your annotation file:

Now you should be able to use the following script:
```
~/RNAseq_assembly_annotation/GAAS/annotation/Tools/bin/maker_gff3manager_JD_v8.pl -f stringtie_with_interpro.gff -b stringtie_annotation.annie -o finalOutputDir
```
That will add the name attribute to the "gene" feature and the description attribute (corresponding to the product information) to the "mRNA" feature into you annotation file. This script may be used for other purpose like to modify the ID value by something more convenient.
The improved annotation is a file named "codingGeneFeatures.gff" inside the finalOutputDir.

## Annotate trinity output

We will have no time to annotate the trinity.fasta today. You can use the same steps that we use previously if you map the transcripts on your genome with for instance the exonerate tool.

Trinity output are often annotated using the [Trinotate pipeline](https://github.com/Trinotate/Trinotate.github.io/wiki). The pipeline is quite similar in what we have done previously but in the trinotate pipeline you are creating proper databases and you do not need to map the transcript before annotation.
The trinotate pipeline was before a set of steps but it has now has being implement

## What's next?

The annoted gff or fasta file can be used in different way in annotation. It can be used to help in the first annotation round to get annotated transcripts to mapped to the genome and then help for the genes annotation. It can also be used after the annotation to complement and improve this one if RNAseq data were not available when the genome annotation was done.
