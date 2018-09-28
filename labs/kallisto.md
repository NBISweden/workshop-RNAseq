---
layout: default
---

# Kallisto and Sleuth 

## Transcript-level quantification with Kallisto

*Kallisto* is an "alignment free" RNA-seq quantification method that runs very fast with a small memory footprint, so that it can be run on most laptops. It is a command-line program that can be downloaded as binary executables for Linux or Mac, or in source code format. For a first insight into the program, read [here](https://liorpachter.wordpress.com/2015/05/10/near-optimal-rna-seq-quantification-with-kallisto/) and for the recently published article, see [here](http://www.nature.com/nbt/journal/vaop/ncurrent/full/nbt.3519.html). There is also a preprint [here](http://arxiv.org/abs/1505.02710).

Kallisto is geared towards quantification on the transcript (isoform) level, rather than the gene level (although the latter can also be done by post-processing Kallisto output.) However, read assignment to transcript isoforms cannot (in general) be done unambiguously, so there is an intrinsic "quantification noise" or variability in this process. Kallisto can thus be run either in a single step (which is very fast) or in "bootstrap" mode (which takes longer, but can be done on several processors in parallel) in order to get uncertainty estimates for the expression levels - a kind of error bars for the quantification process. Running with bootstraps is mandatory if you want to perform differential expression analysis of isoforms with Sleuth (see below). 

Kallisto is primarily meant for quantification of an existing set of FASTA sequences, that is, it does not perform transcript assembly and it cannot quantify the expression of novel transcripts that are not in the transcript index that you provide to it. With that said, you can of course use contigs from an assembly that you have produced in some other program in your Kallisto index. It would also be possible to use the software for e g metagenomics or metatranscriptomics quantification.

## Differential expression with Sleuth

*Sleuth* is a companion package for Kallisto which is used for differential expression analysis of transcript quantifications from Kallisto. While you could use other differential expression packages such as limma or DESeq2 to analyze your Kallisto output, Sleuth also takes into consideration the inherent variability in transcript quantification as explained above. Sleuth also allows the modeling of covariates such as batch, individual, tissue type etc. in the same way as DESeq2/edgeR/limma, which is useful for experimental designs with multiple varying factors. 

Unlike Kallisto, Sleuth is an R package. At the end of a Sleuth analysis, it is possible to view a dynamical graphical presentation of the results where you can explore the differentially expressed transcripts in an intuitive way.

It is still early days for Sleuth and as perhaps mentioned during the lecture on differential expression analysis, it has not been extensively benchmarked against other packages yet. Let's try it on the same A431 data as in the DESeq2 lab!

Since it takes some time to prepare the data, we have pre-computed Kallisto results which you can download from the download area. If you are interested in the steps to get there, please refer to the next section; otherwise jump down to "Running Sleuth" at this point.

## [Only for info, you do not need to do this!] Preparing Sleuth input with Kallisto

Sleuth was designed to work on output from Kallisto (rather than count tables, like DESeq2, or BAM files, like CuffDiff2), so we need to run Kallisto first. (Note that the outputs from other RNA-seq quantifiers like [Salmon](https://github.com/COMBINE-lab/salmon) or [Sailfish](https://github.com/kingsfordgroup/sailfish) can also be used with Sleuth via the new [wasabi](https://github.com/COMBINE-lab/wasabi) package.)

Kallisto is run directly on FASTQ files. We start by downloading the Kallisto software. It can be installed with pip, if you use that (in which case replace all references to "kallisto/kallisto" below with just "kallisto", as the executable will already be in your PATH), but we can also download it:

		wget https://github.com/pachterlab/kallisto/releases/download/v0.42.3/kallisto_mac-v0.42.3.tar.gz

You can of course download a different version if you want to - and if you are using Linux, you need to replace "mac" with "linux" in the command above.

		tar zvxf kallisto_mac-v0.42.3.tar.gz 
		wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-mac64.tar.gz

Again, if you are not using Mac, you need to change the file name above to something appropriate from http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2.

Now we will download and merge human cDNA and ncDNA files from ENSEMBL in order to build a Kallisto transcript index. Note that we can concatenate the two gzipped files without unpacking them first. We use both the protein-coding transcripts and the non-coding ones to be able to capture more of the transcriptome.

	wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
	wget ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
	cat Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.ncrna.fa.gz > Homo_sapiens.GRCh38.rna.fa.gz

Now we can build the transcriptome index. Let's also time it to get a sense of how long it takes.
	
	time kallisto/kallisto index -i hsGRCh38_kallisto Homo_sapiens.GRCh38.rna.fa.gz

It should take less than 10 minutes.

Next, copy the FASTA files from uppmax ``/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/diffExp/FASTQ/`` .
When that is done, it's time for quantifying the FASTQ files against our Kallisto index with bootstrapping for later use in Sleuth. You could do that one by one, with a command like

    time kallisto/kallisto quant -i hsGRCh38_kallisto -t 4 -b 100 7_111116_AD0341ACXX_137_1_index1_1.fastq.gz 7_111116_AD0341ACXX_137_1_index1_2.fastq.gz -o sample1
	
or in a bash loop:

    for i in {1..12}; do time kallisto/kallisto quant -i hsGRCh38_kallisto -t 4 -b 100 7_111116_AD0341ACXX_137_${i}_index${i}_1.fastq.gz 7_111116_AD0341ACXX_137_${i}_index${i}_2.fastq.gz -o sample${i}; done

In this example, we put "-t 4" so we can use up to four processors in the bootstrapping. You may want to modify this value according to the machine you are working on. If you wanted to run Kallisto without bootstraps and just get expression values on a pair of FASTQ files, you would run something like

    kallisto/kallisto quant -i hsGRCh38_kallisto <FILE1>.fastq <FILE2>.fastq -o <OUTPUT_DIR_NAME>

Running Kallisto on all the 12 samples with 100 bootstraps may take an hour or so, depending on your machine and settings. The time on a MacBook Pro with four threads and 16 Gb of RAM was ... minutes.

## Running Sleuth

Here we give an example workflow for a DE analysis in Sleuth based on the A431 data that we are using for all the DE analysis labs. Start by copy the results from uppmax ``/proj/uppstore2017171/courses/RNAseqWorkshop/downloads/diffExp/kallisto_results.tar.gz``. Download and extract the whole folder and make a note of where it is.

The Sleuth analysis is done entirely in R, so start your R environment and begin by installing the dependencies. This only needs to be done the first time, of course.

		source("http://bioconductor.org/biocLite.R")
		biocLite("rhdf5")
		install.packages("devtools") 
		library("devtools")
		devtools::install_github("pachterlab/sleuth")
 
Now load the package and use a function that we borrowed from the Sleuth documentation for connecting ENSEMBL transcript names to common gene names. This will turn out to be useful at the end, when we look at the dynamic visualization of the results.

		library("sleuth")

		tx2gene <- function(){
		mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
		t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                	"external_gene_name"), mart = mart)
		t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     	ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
		return(t2g)
		}

		t2g <- tx2gene()

The actual Sleuth analysis starts by defining the path to the directory where your Kallisto output folders are. You will need to replace "/PATH/TO/YOUR/FOLDER" with the actual path on your machine. 

		base_dir <- "PATH/TO/YOUR/FOLDER"
		
Then we tell the program about the samples that we have. We stored the Kallisto results for the 12 samples in directories called "sample1", "sample2", ..., "sample12".

		samples <- paste0("sample", 1:12)
		kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

Now it's time to fill in metadata about the samples. We can use a similar assignment as in the DESeq2 exercise:

    s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = rep(c("ctrl", "t2h", "t6h", "t24h"), each=3), stringsAsFactors=FALSE)

Again, if there were other experimental factors involved, these could have been modelled here as well. If you want to look at such an example, you might want to refer to the [beta version of this exercise](http://scilifelab.github.io/courses/rnaseq/labs/kallisto) that was given in October 2015. In that version, we did not use A431 data but rather a prostate cancer data set where the two experimental factors were (1) the individual that the sample came from, (2) tumor or normal tissue.
		
Back to the present data! The next command will read the Kallisto output files, connect them with metadata, and set up a linear model for analyzing the expression data.
 
		so <- sleuth_prep(s2c, ~timepoint, target_mapping = t2g)

Next we fit the linear model and test for one of the model coefficients. In this case we test the 24h time point versus the control.

		so <- sleuth_fit(so)
		so <- sleuth_wt(so, which_beta="timepointt24h") 

Now we should be able to visualize the results:

		sleuth_live(so)
	
There are lots of things to look at here - explore according to your interests! Some things you might try are e.g. the PCA and sample heatmap options in the map menu, the test table in the analyses menu (which contains a ranked list of the differentially expressed genes), or the gene view in the same menu.

If you want to delve further into time series analysis with Sleuth (after all, we have just compared two time points here, whereas we have four in all), you might want to read this [excellent blog post](http://nxn.se/post/134227694720/timecourse-analysis-with-sleuth) by Valentine Svensson. Note that Sleuth is still under development, so some of the commands may be a bit different.

