---
layout: default
title:  'R-intro'
---

# Very quick introduction to R

In several exercises in this course, we use the statistical
programming environment R. Here follows a very quick introduction,
primarily for those who have never used R before. You can go through
this now, or return to it when you first encounter R in the exercises.

First start R on your computer. How to do this depends on your
operating system. If you are on a Linux or Mac OS X system, you would
typically execute the command ``R`` at the shell prompt.

In this course, we are running R on a Linux server. If you would like
to install R on your own computer, you can find the appropriate
download for your system [here](http://ftp.acc.umu.se/mirror/CRAN/).

When R starts, you should see a text message stating the version of
R you are running and some further information, followed by a
command-line “prompt” ( ``>``, a greater-than sign). The prompt means that
R is waiting for you to type a command. To get a feel for how this
works, try out some arithmetic::

	 > 2 + 3
	 [1] 5
	 > 5 * 4 + 10
	 [1] 30

And some function calls

	> abs(-5)
	[1] 5
	> sum(1,5,10)
	[1] 16

To find out how to use a function, type its name preceded by a question mark::

	> ? sin

This will bring up some help documentation for the function. You can
use the arrow keys to scroll the help text up and down. Press q to get
back to the R prompt.

Now try this:
	
	> a <- 10

This command created an object called a. Objects are an important
concept in R (as in many other programming languages), and we will be
creating more of them in the RNA-seq exercises. We can inspect an
object by just typing its name:
	
	> a
	[1] 10

We can also change the value of an object that we’ve created:

	> a <- 2 * a	
	> a
	[1] 20

The object *a* created above is a vector with a single element. To
create a vector with several elements, you can use the function *c*:


	> b <- c(1, 2, 10)
	> b
	[1]  1  2 10

Or the colon operator::

	> 1:10
	[1]  1  2  3  4  5  6  7  8  9 10

A matrix can be created with the function *cbind*:

	> b <- cbind(1:10, 101:110)
	> b
       [,1] [,2]
	[1,]    1  101
	[2,]    2  102
	[3,]    3  103
	[4,]    4  104
	[5,]    5  105
	[6,]    6  106
	[7,]    7  107
	[8,]    8  108
	[9,]    9  109
	[10,]   10  110

We can then use indices to access selected elements of the matrix:

	> b[1,]
	[1]   1 101
	> b[, 2]
	[1] 101 102 103 104 105 106 107 108 109 110
	> b[c(5,8), 2]
	[1] 105 108

To end your R session, use the functions *quit* or *q*:
 
	> q()
	Save workspace image? [y/n/c]: n

If you wish, you can save a workspace image. This means that the state of your R session is written to a file, which will be loaded the next time
you start R in the same working directory. This is good if you want to continue the session at a later time. If not, just answer *n* (for no).

You can find manuals for R and more information on the [R web site](http://www.r-project.org/).

## Working with R packages

R contains many functions, in particular for statistics and
plotting. There are numerous additional functions and data sets,
contained in R packages that are not installed by default.

Many such extension packages can be obtained from the
[Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/web/packages/index.html).
We will also download packages from [BioConductor](http://www.bioconductor.org), an alternative
repository focused on R packages for bioinformatics.

### Installing packages from CRAN

To install a CRAN package in R, use the install.packages() function. This simple command downloads the package from a specified repository (by default, CRAN) and installs it on your machine. For example, we can install a useful package called *gplots*:

	> install.packages("gplots")

After installation, the package is stored in a location on your computer called the R package library. Packages are loaded from this library using the function *library*:

	> library("gplots")

One useful function in the *gplots* package is *heatmap.2*, which plots heatmaps with dendrograms. Since we have installed and loaded the *gplots* package, we can now use this function. Let's try it with a matrix of random, normally distributed values:

	> heatmap.2( matrix(rnorm(180), ncol=6) )

Note that you only need to install a package once. The next time you start R, you can run *library("gplots")* directly, because the package has already been installed. If you upgrade to a new version of R, you may have to reinstall packages.

### Installing packages from BioConductor

To install packages from BioConductor, you first have to load a function *biocLite.R* from BioConductor and then use that function to install the packages.
That is done in two lines of code:

	> source("https://bioconductor.org/biocLite.R")
	> biocLite("Name of package")

For example, let's install the package *DESeq2*, which provides functions for RNA-seq differential expression analysis:
	
	> source("https://bioconductor.org/biocLite.R")
	> biocLite("DESeq2")

You should now be able to load the package *DESeq2*:

	> library("DESeq2")

### Installing packages from GitHub

Sometimes you want to install packages that are still under development. Those are mostly found on [GitHub](https://github.com).

The example below shows how you can switch between a stable version and a development version of a package, in this case the package *ggplot2*.
This example is provided for reference and you don't have to understand it (or try it) for completing the RNA-seq labs.

First, we install and load the regular (stable) version of the *ggplot2* package from CRAN:

	> install.packages("ggplot2")
	> library(ggplot2)

To see which packages are loaded, and their versions, use the function *sessionInfo*:

	> sessionInfo()

Note the version number of *ggplot2* reported by *sessionInfo*.

There is no entirely safe way to unload packages in R, so to work with a different package version, it is safest to quit and restart R. To quit:

	> q()
	Save workspace image? [y/n/c]: n

After starting R again, we install and load the *devtools* package, which provides functions for working with development versions of packages:

	> install.packages("devtools")
	> library(devtools)

Next we install the development version of *ggplot2* from GitHub:

	> dev_mode(on=TRUE)
	> install_github("hadley/ggplot2")

Now we can load and work with the the development version of *ggplot2*:

	> library(ggplot2)
	> sessionInfo()

Check the output from *sessionInfo* to verify that you now get a different version number for *ggplot2*!

When finished working with development versions, do:

	> dev_mode(on=FALSE)

Note that if you only want to install from GitHub, you can ignore the *dev_mode* function calls above.
