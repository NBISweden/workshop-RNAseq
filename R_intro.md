---
layout: default
title:  'R instructions'
---

# R instructions

You need to have R installed in order to do the exercises. It will also be useful to have Rstudio installed, to run all your commands. Both of these works on computers running Linux, Windows and Macintosh operative systems. R-studio is a set of tools as well as an editor that facilitates the use of R and in many ways create a superior environment to integrate directly with R. Over the last years it has become a very popular tool and in many ways become a de-facto standard for working with R.

Note that on same operative systems it will be easier to install and run R and R-studio if you are administrator of your own computer and hence are allowed to install software on your machine. If you do not have these privileges please ask your system administrator to install the latest version of R and R-studio.

### Install R
1. Go to [CRAN](https://cran.rstudio.com/)
2. Click on the link corresponding to your operative system
3. Download the recommended files for your system.
4. Run the installer or move the downloaded files to suitable place on your computer.

##### Details for Windows  
Once you clicked on the “Download R for Windows” you will come to a new web page where you will have a set of options. Click on the first link named “base”. This will take you to the web page entitled “R-3.3.3 for Windows” (the numbers might be different depending on what is the most recent version) where you can download the “R-3.3.3-win.exe” that will can be run to install R on your computer.

##### Details for Macintosh

Once you clicked on the “Download R for Macintosh” you will come to a new web page where you will have a set of options. Unless you have an old version of your operative system you should select the first link named “R-3.3.3.pkg”  (the numbers might be different depending on what is the most recent version) that will download R to your computer. If you are not sure what version you are running click the apple on the top left of your screen and select “About this mac” (Om den här datorn). 

You can then double-click the downloaded package that will prompt you with some questions for installation details. Stick with the default settings and you should be fine.

##### Details for Linux

Once you clicked on the “Download R for Linux” you will come to a new web page where you can select the linux version you use. On most distributions this will be via a software install system like yum or apt-get. If you run this make sure that you update your information to the installer first, otherwise you might end up installing at outdated version of R. For some systems you might need to install not only r-base, but also r-devel or you will lack important features of your R installation.

### Install R-studio

Go to the web page [rstudio](https://www.rstudio.com/products/rstudio/download/) download the installer corresponding to your operative system. Unpack the installer and install the app on a suitable place on your system.

### Useful R online resources <a name="R_tutorials"></a>

If you need to touch up on your R-skills before the course start, the following links goes to some useful material on using R.

* [Best first R tutorial](https://www.nceas.ucsb.edu/files/scicomp/Dloads/RProgramming/BestFirstRTutorial.pdf)
A nice self learn tutorial to R, introducing many central concepts to R.
* [A short introduction to R](https://cran.r-project.org/doc/contrib/Torfs+Brauer-Short-R-Intro.pdf)
A very short intro to using R.
* [An introduction to R](https://cran.r-project.org/doc/manuals/r-release/R-intro.html)  
A fairly comprehensive document on R. As a beginner one can start with Appendix A that is a short practical session.
* [The art of R programming](http://heather.cs.ucdavis.edu/~matloff/132/NSPpart.pdf)
A pdf copy of a book that deals with R as a programming language. Is a great source of information for anyone that wants to use R not only as a statistical analysis tools, but also use it as a more general programming language.
* [R for data science](http://r4ds.had.co.nz/)
The basics of getting data into R, clean the data do your analysis. One of the authors of this book Hadley Wickham is also behind R-studio and has published several books on different aspects of using R, many of them available for free via his web page. He has also created many packages for R that facilitate structured data analysis. His most popular packages is ggplot2 for creating beautiful graphics with limited set of commands. Besides ggplot2 he has created many other packages that make use of novel objects and often use slightly different R syntax to interact with them. This means that code using them will look different to most other R code. If you like this way of working, there is a whole set of packages often referred to as the [tidyverse](https://blog.rstudio.org/2016/09/15/tidyverse-1-0-0) that offers a more unified interaction with objects and functions under this approach to data analysis.

