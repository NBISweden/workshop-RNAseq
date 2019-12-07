---
layout: default
title:  'Linux Intro'
---

# Linux Introduction
**NOTE:** in syntax examples, the dollar sign ($) is not to be printed.
The dollar sign is usually an indicator that the text following it should be typed in a terminal window.

## 1. Connecting to UPPMAX
The first step of this lab is to open a ssh connection to UPPMAX.
You will need a ssh program to do this:

On Linux: it is included by default, named **Terminal**.

On OSX: it is included by default, named **Terminal**.

On Windows: [Google MobaXterm](http://bit.ly/19yaQOM) and download it.

Fire up the available ssh program and enter the following (replace **username** with your uppmax user name).
-Y means that X-forwarding is activated on the connection, which means graphical data can be transmitted if a program requests it, i.e. programs can use a graphical user interface (GUI) if they want to.

```
$ ssh -Y username@rackham.uppmax.uu.se
```

and give your password when prompted.
As you type, nothing will show on screen.
No stars, no dots.
It is supposed to be that way.
Just type the password and press enter, it will be fine.

Now your screen should look something like this:
```
dahlo@dahlo-xps ~ $ ssh -Y dahlo@rackham.uppmax.uu.se
Last login: Fri May 18 15:03:59 2018 from micro046.icm.uu.se
 _   _ ____  ____  __  __    _    __  __
| | | |  _ \|  _ \|  \/  |  / \   \ \/ /   | System:    rackham4
| | | | |_) | |_) | |\/| | / _ \   \  /    | User:      dahlo
| |_| |  __/|  __/| |  | |/ ___ \  /  \    | 
 \___/|_|   |_|   |_|  |_/_/   \_\/_/\_\   | 

###############################################################################

        User Guides: http://www.uppmax.uu.se/support/user-guides
        FAQ: http://www.uppmax.uu.se/support/faq

        Write to support@uppmax.uu.se, if you have questions or comments.


dahlo@rackham4 ~ $ 
```
## 2. Getting a node of your own

Usually you would do most of the work in this lab directly on one of the login nodes at uppmax, but we have arranged for you to have one core each to avoid disturbances.
This was covered briefly in the lecture notes.

```
$ salloc -A snic2019-8-3 -t 07:00:00 -p core -n 1 --no-shell --reservation=snic2019-8-3_4 &
```

check which node you got (replace **username** with your uppmax user name)

```
$ squeue -u username
```

should look something like this

```
dahlo@rackham2 work $ squeue -u dahlo
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
           3132376      core       sh    dahlo  R       0:04      1 r292
dahlo@rackham2 work $ 
```

where **r292** is the name of the node I got (yours will probably be different).
Note the numbers in the Time column.
They show for how long the job has been running.
When it reaches the time limit you requested (7 hours in this case) the session will shut down, and you will lose all unsaved data.
Connect to this node from within uppmax.

```
$ ssh -Y r292 
```

**Note:** there is a uppmax specific tool called jobinfo that supplies the same kind of information as squeue that you can use as well (```$ jobinfo -u username```).

## 3. Moving and Looking Around
It is good to know how to move around in the file system.
I'm sure you all have experienced this using a graphical user interface (**GUI**) before, Windows Explorer in Windows and Finder in OSX.
Using the command line can be confusing at first, but the more your do it, the easier it gets.

When you connect to UPPMAX, you will start out in your home folder.
The absolute path to your home folder is usually **/home/\<username\>**

Start with looking at what you have in your home folder.
The command for this ls, and it stand for **L**i**S**t (list).

```
$ ls -l
```

You should see something similar to:

```
[12:56:32] dahlo@rackham1 ~ $ ls -l
total 4384
drwxr-xr-x  3 dahlo b2014068    2048 Apr 27  2017 archive
-rw-rw-r--  1 dahlo dahlo      49042 Sep 23  2016 bad.png
-rw-r--r--  1 dahlo root     1704126 Feb 16  2015 CORRUPTED_FILE_LIST.20150216
drwxr-xr-x  2 dahlo dahlo       2048 Mar 18  2016 data
drwxr-x---  3 dahlo dahlo       2048 Feb 15  2017 glob
drwxrwxr-x  3 dahlo dahlo       2048 Feb 10  2015 glob2
-rw-rw-r--  1 dahlo dahlo      60944 Sep 23  2016 good.png
drwxr-xr-x  4 dahlo b2014209    2048 Oct 30  2014 igv
drwxrwxr-x  5 dahlo dahlo       2048 Sep 20  2016 ngsintro
drwxr-xr-x  5 dahlo dahlo       2048 Jun  3  2014 perl
drwxrwxr-x  5 dahlo dahlo       2048 Nov 22  2012 perl5
drwx--S---  2 dahlo dahlo       2048 May  4  2010 private
drwxr-xr-x  6 dahlo dahlo       2048 May 17  2016 R
lrwxrwxrwx  1 dahlo dahlo         41 Nov 14  2013 thindrives -> /var/opt/thinlinc/sessions/dahlo/6/drives
drwxr-xr-x 26 dahlo dahlo       4096 May 18 10:43 uppmaxScripts
drwxrwxr-x  5 dahlo dahlo    2201600 May 14 14:02 work
[12:57:36] dahlo@rackham1 ~ $ 
```

We can start off with moving to the course's **proj folder**.
Each project has its own **proj folder** and **nobackup folder** that contains all the files that belong to the project.
The proj folder is backed up, so if a file is deleted by mistake it can be retrieved from the backups, and this is where only compressed raw data and scripts are usually stored.
The nobackup folder is as the name implies not backed up, and this is where all the other files belonging to the project is kept, like everything that is created by programs that you run.
Backup is expensive, so we don't want to waste it on temporary files generated by your analysis programs, which could be re-run.
The nobackup folder is the folder you will keep most of your files in.

As seen in the lecture, the command for moving around is **cd**.
The command stands for **C**hange **D**irectory and does exactly that.
It is the equivalent of double clicking a folder in a GUI.
To enter the course project's proj folder, simply type

```
$ cd /proj/snic2019-8-3
```

Have a look around and see which folders there are.
This is the folder where all the unreplacable files belonging to a project would go, like the raw sequencing data that is delivered by the sequencing platforms.
Since we'll be creating a bunch of temporary files during these labs, we will have to do that in the nobackup folder.
Go there using the cd command.

```
$ cd nobackup
```

We can easily see that this is a relative path, since it does not start with a **/** That means that this command will only work when you are standing in your project folder.
If you are standing somewhere else and say that you want to enter a folder named nobackup, the computer will tell you that there is no folder named nobackup where you are located at the moment.
The absolute path to the course's nobackup folder would be **/proj/snic2019-8-3/nobackup**

It is the exact same things as if you are using a GUI.
If you are standing on your desktop, you can double click a folder which is located on your desktop.
But if you are standing in another folder, you can't double click on the same folder, because it is just not located in that folder.
You have to move to your desktop first, and then double click it.

If you look around in the nobackup folder, you probably only have a folder called **private**.
This is how every nobackup folder looks before you start putting files there.

Typing ls -l all the time is.. more annoying than one would think, so someone came up with the good idea to add a shortcut here.
If you type **ll**, it is the same as typing **ls -l** Use it from now on.

Now we have practiced moving around and looking at what we have in folders.
The next step will show you how to do the same thing, but without the moving around part.
If we want to look at what we have in our home folder, while standing in the course's project folder, we type ```$ ll /home/<username>/``` and remember to substitute **\<username\>** with your own user name.

```
$ ll /home/valent/
```

Since most programmers are lazy (efficient), there is a shortcut to your home folder so that you don't have to write it all the time.
If you write **~/** it means the same as if you would write **/home/\<username\>/**
Try using it with ls:

```
$ ll /home/<username>
```

or

```
$ ll ~/
```

## 4. Copying files needed for laboratory
To be able to do parts of this lab, you will need some files.
To avoid all the course participants editing the same file all at once, undoing each other's edits, each participant will get their own copy of the needed files.

The files are located in the folder **/sw/share/compstore/courses/ngsintro/linux_tutorial**

or they can be downloaded if you are not on uppmax at the moment, [files.tar.gz](files/linux-intro/files.tar.gz) (instruction on how to download further down)

For structures sake, first create a folder named as your username in the course's nobackup folder, and a folder called **linux_tutorial** inside that folder, where you can put all your lab files.

This can be done in 2 ways:

```
$ mkdir /proj/snic2019-8-3/nobackup/<username>
$ mkdir /proj/snic2019-8-3/nobackup/<username>/linux_tutorial
```

or

```
$ mkdir -p /proj/snic2019-8-3/nobackup/<username>/linux_tutorial
```

The reason for this is that Linux will not like it if you try to create the folder linux_tutorial inside a folder (the one named like your username) that does not exist yet.
Then you have the choice to either first create the one named like your username (the first way), or to tell Linux to create it for you by giving it the **-p** option (the second way).

Next, copy the lab files to this folder.
**-r** means recursively, which means all the files including sub-folders of the source folder.
Without it, only files directly in the source folder would be copied, NOT sub-folders and files in sub-folders.

**NOTE:** Remember to tab-complete to avoid typos and too much writing.

```
Ex: cp -r <source folder> <destination folder>

$ cp -r /sw/share/compstore/courses/ngsintro/linux_tutorial/* /proj/snic2019-8-3/nobackup/<username>/linux_tutorial
```

If you are unable to copy the files on uppmax, you can download the files instead of copying them.
This is done with the command **wget** (web get).
It works kind of the same way as the cp command, but you give it a source URL instead of a source file, and you specify the destination by giving it a prefix, a path that will be appended in front on the file name when it's downloaded.
I.e. if you want to download the file http://somewhere.com/my.file and you give it the prefix ~/analysis/, the downloaded file will be saved as ~/analysis/my.file

```
Ex: wget -P <destination prefix> <source URL>

$ wget -P /proj/snic2019-8-3/nobackup/<username>/linux_tutorial {{site.url}}/ngsintro/common/files/linux_additional-files.zip
```

## 5. Unpack Files
Go to the folder you just copied and see what is in it.

**NOTE: Remember to tab-complete to avoid typos and too much writing.**

```
$ cd /proj/snic2019-8-3/nobackup/<username>/linux_tutorial

$ ll
```

tar.gz is a file ending give to compressed files, something you will encounter quite often.
Compression decreases the size of the files which is good when downloading, and it can take thousands of files and compress them all into a single compressed file.
This is both convenient for the person downloading and speeds up the transfer more than you would think.

To unpack the **files.tar.gz** file use the following line while standing in the newly copied linux_tutorial folder.

```
$ tar -xzvf files.tar.gz
```

The command will always be the same for all tar.gz files you want to unpack.
-xzvf means e<b>X</b>tract from a **Z**ipped file, **V**erbose (prints the name of the file being unpacked), from the specified **F**ile (f must always be the last of the letters).

Look in the folder again and see what we just unpacked:
```
[valent@milou2 linux_tutorial]$ ll -la
total 512
drwxrwsr-x 12 valent snic2019-8-3   2048 Sep 24 13:19 .
drwxrwsr-x  6 valent snic2019-8-3   2048 Sep 24 13:19 ..
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 a_strange_name
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 backed_up_proj_folder
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 external_hdd
-rwxrwxr-x  1 valent snic2019-8-3  17198 Sep 24 13:19 files.tar.gz
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 important_results
drwxrwsr-x  2 valent snic2019-8-3 129024 Sep 19  2012 many_files
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 old_project
-rwxrwxr-x  1 valent snic2019-8-3      0 Sep 24 13:19 other_file.old
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 part_1
drwxrwsr-x  2 valent snic2019-8-3   2048 Sep 19  2012 part_2
drwxrwsr-x  2 valent snic2019-8-3   2048 Jan 28  2012 this_has_a_file
drwxrwsr-x  2 valent snic2019-8-3   2048 Jan 28  2012 this_is_empty
-rwxrwxr-x  1 valent snic2019-8-3      0 Sep 19  2012 useless_file
[valent@milou2 linux_tutorial]$ 
```
## 6. Copying and Moving Files
Let's move some files.
Moving files might be one of the more common things you do, after cd and ls.
You might want to organize your files in a better way, or move important result files to the project folder, who knows?

We will start with moving our important result to a backed up folder.
When months of analysis is done, the last thing you want is to lose your files.
Typically this would mean that you move the final results to your project folder.

In this example we want to move the result file only, located in the folder `important_results`, to our fake project folder, called `backed_up_proj_folder`.
The syntax for the move command is:

```
$ mv <source> <destination>
```

First, take a look inside the important_results folder:

```
[valent@milou2 linux_tutorial]$ ll important_results/
total 0
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 dna_data_analysis_result_file_that_is_important-you_should_really_use_tab_completion_for_file_names.bam
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 temp_file-1
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 temp_file-2
[valent@milou2 linux_tutorial]$
```
You see that there are some unimportant temporary files that you have no interest in.
Just to demonstrate the move command, I will show you how to move one of these temporary files to your backed up project folder:

```
$ mv important_results/temp_file-1 backed_up_proj_folder/
```

**Now do the same, but move the important DNA data file!**

Look in the backed up project folder to make sure you moved the file correctly.

```
[valent@milou2 linux_tutorial]$ ll backed_up_proj_folder/
total 0
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 dna_data_analysis_result_file_that_is_important-you_should_really_use_tab_completion_for_file_names.bam
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 last_years_data
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 temp_file-1
[valent@milou2 linux_tutorial]$
```
Another use for the move command is to rename things.
When you think of it, renaming is just a special case of moving.
You move the file to a location and give the file a new name in the process.
The location you move the file **to** can of very well be the same folder the file already is in.
To give this a try, we will rename the folder `a_strange_name` to a better name.

```
$ mv a_strange_name a_better_name
```

Look around to see that the name change worked.

```
[valent@milou2 linux_tutorial]$ mv a_strange_name a_better_name
[valent@milou2 linux_tutorial]$ ll
total 448
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 a_better_name
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 24 13:40 backed_up_proj_folder
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 external_hdd
-rwxrwxr-x 1 valent snic2019-8-3  17198 Sep 24 13:36 files.tar.gz
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 24 13:40 important_results
drwxrwsr-x 2 valent snic2019-8-3 129024 Sep 19  2012 many_files
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 old_project
-rwxrwxr-x 1 valent snic2019-8-3      0 Sep 24 13:36 other_file.old
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 part_1
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 part_2
drwxrwsr-x 2 valent snic2019-8-3   2048 Jan 28  2012 this_has_a_file
drwxrwsr-x 2 valent snic2019-8-3   2048 Jan 28  2012 this_is_empty
-rwxrwxr-x 1 valent snic2019-8-3      0 Sep 19  2012 useless_file
[valent@milou2 linux_tutorial]$
```

Sometimes you don't want to move things, you want to copy them.
Moving a file will remove the original file, whereas copying the file will leave the original untouched.
An example when you want to do this could be that you want to give a copy of a file to a friend.
Imagine that you have a external hard drive that you want to place the file on.
The file you want to give to your friend is data from last years project, which is located in your backed_up_project_folder, `backed_up_proj_folder/last_years_data`

As with the move command, the syntax is

```
$ cp <source> <destination>

$ cp backed_up_proj_folder/last_years_data external_hdd/
```

Take a look in the external_hdd to make sure the file got copied.
```
[valent@milou2 linux_tutorial]$ cp backed_up_proj_folder/last_years_data external_hdd/
[valent@milou2 linux_tutorial]$ ll external_hdd/
total 0
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 24 13:46 last_years_data
[valent@milou2 linux_tutorial]$
```
## 7. Deleting Files
Sometimes you will delete files.
Usually this is when you know that the file or files are useless to you, and they only take up space on your hard drive or UPPMAX account.
To delete a file, we use the **R**e<b>M</b>ove command, **rm**.
Syntax:

```
$ rm <file to remove>
```

If you want, you can also specify multiple files a once, as many as you want!

```
$ rm <file to remove> <file to remove> <file to remove> <file to remove> <file to remove>
```

**IMPORTANT: There is no trash bin in Linux.
If you delete a file, it is gone.
So be careful when deleting stuff..**

Try it out by deleting the useless file in the folder you are standing in.
First, look around in the folder to see the file.

```
[valent@milou2 linux_tutorial]$ ll
total 448
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 a_better_name
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 24 13:40 backed_up_proj_folder
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 24 13:46 external_hdd
-rwxrwxr-x 1 valent snic2019-8-3  17198 Sep 24 13:36 files.tar.gz
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 24 13:40 important_results
drwxrwsr-x 2 valent snic2019-8-3 129024 Sep 19  2012 many_files
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 old_project
-rwxrwxr-x 1 valent snic2019-8-3      0 Sep 24 13:36 other_file.old
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 part_1
drwxrwsr-x 2 valent snic2019-8-3   2048 Sep 19  2012 part_2
drwxrwsr-x 2 valent snic2019-8-3   2048 Jan 28  2012 this_has_a_file
drwxrwsr-x 2 valent snic2019-8-3   2048 Jan 28  2012 this_is_empty
-rwxrwxr-x 1 valent snic2019-8-3      0 Sep 19  2012 useless_file
[valent@milou2 linux_tutorial]$
```
Now remove it:

```
$ rm useless_file
```

Its not only files you can remove.
Folders can be removed too.
There is even a special command for removing folders, **rmdir**.
They work similar to **rm**, except that they can't remove files.
There are two folders, `this_is_empty` and `this_has_a_file`, that we now will delete.

```
$ rmdir this_is_empty

$ rmdir this_has_a_file
```

If you look inside `this_has_a_file`:

```
[valent@milou2 linux_tutorial]$ ll this_has_a_file
total 0
-rwxrwxr-x 1 valent snic2019-8-3 0 Jan 28  2012 file
[valent@milou2 linux_tutorial]$
```

you will see there is a file in there! Only directories that are completely empty can be deleted using rmdir.
To be able to delete `this_has_a_file`, either delete the file manually and then remove the folder

```
$ rm this_has_a_file/file

$ rmdir this_has_a_file
```

or delete the directory recursively, which will remove `this_has_a_file` and everything inside:

```
$ rm -r this_has_a_file
```

## 8. Open files
So what happens if you give your files bad names like 'file1' or 'results'? You take a break in a project and return to it 4 months later, and all those short names you gave your files doesn't tell you at all what the files actually contain.
Of course, this should never happened, because you should **ALWAYS** name your files so that you definetly know what they contain.
But lets say it did happen.
Then the only way out is to look at the contents of the files and try to figure out if it is the file you are looking for.

Now, we are looking for that really good script we wrote a couple of months ago in that other project.
Look in the project's folder, **old_project**, and find the script.

```
[valent@milou2 linux_tutorial]$ ll old_project/
total 96
-rwxrwxr-x 1 valent snic2019-8-3 39904 Sep 19  2012 a
-rwxrwxr-x 1 valent snic2019-8-3     0 Sep 19  2012 stuff_1
-rwxrwxr-x 1 valent snic2019-8-3  1008 Sep 19  2012 the_best
[valent@milou2 linux_tutorial]$
```

Not so easy with those names..
We will have to use less to look at the files and figure out which is which.
Syntax for less:

```
$ less <filename>
```

Have a look at **the_best**, that must be our script, right?

```
$ less old_project/the_best
```

**(press q to close it down, use arrows to scroll up/down)**

I guess not.
Carrot cakes might be the bomb, but they won't solve bioinformatic problems.
Have a look at the file **a** instead.

That's more like it!

Now imagine that you had 100s of files with weird names, and you really needed to find it..
Lesson learned: name your files so that you know what they are! And don't be afraid to create folders to organize files.

Another thing to think about when opening files in Linux is which program should you open the file in? The programs we covered during the lectures are **nano** and **less**.
The main difference between these programs in that **less can't edit files**, only view them.
Another difference is that **less doesn't load the whole file** into the RAM memory when opening it.
So, why care about how the program works? I'll show you why.
This time we will be opening a larger file, located in the course's project folder.
It's 65 megabytes, so it is a tiny file compared with bio-data.
Normal sequencing files can easily be 100-1000 times larger than this.

First, open the file with nano.
Syntax:

```
$ nano <filename>

$ nano /sw/share/compstore/courses/ngsintro/linux_additional-files/large_file
```

**(press ctrl+x to close it down, user arrows to scroll up/down)**

Is the file loaded yet? Now take that waiting time and multiply it with 100-1000.
Now open the file with less.
Notice the difference?

**Head** and **tail** works the same was as less in this regard.
They don't load the whole file into RAM, they just take what they need.

To view the first rows of the large file, use head.
Syntax:

```
$ head <filename>

$ head /sw/share/compstore/courses/ngsintro/linux_additional-files/large_file
```

Remember how to view an arbitrary number of first rows in a file?

```
$ head -n <number of rows to view> <filename>

$ head -n 23 /sw/share/compstore/courses/ngsintro/linux_additional-files/large_file
```

The same syntax for viewing the last rows of a file with tail:

```
$ tail <filename>

$ tail /sw/share/compstore/courses/ngsintro/linux_additional-files/large_file

$ tail -n <number of rows to view> <filename>

$ tail -n 23 /sw/share/compstore/courses/ngsintro/linux_additional-files/large_file
```

## 9. Wildcards
Sometimes (most of the time really) you have many files.
So many that it would take you a day just to type all their names.
This is where **wildcards** saves the day.
The wildcard symbol in Linux is the star sign, `*` , and it means literally **anything**.
Say that you want to move all the files which has names starting with `sample_1_` and the rest of the name doesn't matter.
You want all the files belonging to `sample_1`.
Then you could use the wildcard to represent the rest of the name:

**(don't run this command, it's just an example)**

```
$ mv  sample_1_*  my_other_folder
```

We can try it out on the example files I have prepared.
There are two folder called **part_1** and **part_2**.
We want to collect all the **.txt** files from both these folders in one of the folders.
Look around in both the folders to see what they contain.

```
[valent@milou2 linux_tutorial]$ ll part_1/
total 0
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 file_1.txt
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 file_2.txt
[valent@milou2 linux_tutorial]$ ll part_2
total 0
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 file_3.txt
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 file_4.txt
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 garbage.tmp
-rwxrwxr-x 1 valent snic2019-8-3 0 Sep 19  2012 incomplete_datasets.dat
[valent@milou2 linux_tutorial]$
```

We see that **part_1** only contains .txt files, and that **part_2** contains some other files as well.
The best option seem to be to move all .txt files from **part_2** info **part_1**.

```
$ mv part_2/*.txt part_1/
```

The wildcard works with most, if not all, Linux commands.
We can try using wildcards with **ls**.
Look in the folder `many_files`.
Yes, there are ~1000 .docx files in there.
But not only .docx files..
There are a couple of **.txt** files in there as well.
Find out which numbers they have.
Try to figure out the solution on your own.

<details>
  <summary>I give up, tell me how to do it!</summary><p>
 
```
$ ll many_files/*.txt
```
</p></details>

## 10. Utility Commands
Ok, the last 2 commands now.

**top** and **man**.

Top can be useful when you want to look at which programs are being run on the computer, and how hard the computer is working.
Type **top** and have a look.

```
$ top
```

![](files/linux-intro/top.jpg)

**(press q to close it down)**

Each row in top corresponds to one program running on the computer, and the column describe various information about the program.
**The right-most column** shows you which program the row is about.

There are mainly 2 things that are interesting when looking in top.
The first is **how much cpu** each program is using.
I have marked it with **blue** in the picture.
If you are doing calculations, which is what bioinformatics is mostly about, the cpu usage should be high.
The numbers in the column is how many percent of a core the program is running.
If you have a computer with 8 cores, like the UPPMAX computers, you can have 8 programs using 100% of a core each running at the same time without anything slowing down.
As soon as you start a 9th program, it will have to share a core with another program and those 2 programs will run at half-speed since a core can only work that fast.
In the image, one instance of the program ssh is using 78.4% of a core.

The areas marked with **red** is describing **how much memory** is being used.
The area in the top describes the overall memory usage.
Total tells you how much memory the computer has, used tells you how much of the memory is being used at the moment, and free tells you how much memory is free at the moment.
Total = Used + Free

The **red column** tells you how much memory each program uses.
The numbers mean how many percent of the total memory a program uses.
In the image, the program celeryd is using 0.3% of the total memory.

A warning sign you can look for in top is when you are running an analysis which seems to take forever to complete, and you see that there is almost no cpu usage on the computer.
That means that the computer is not doing any calculation, which could be bad.
If you look at the memory usage at the same time, and see that it's maxed out (used â‰ƒ total), you can more or less abort the analysis.

When the memory runs out, the computer more or less stops.
Since it can't fit everything into the RAM memory, it will start using the hard drive to store the things it can't fit in the RAM.
Since the hard drive is ~1000 times slower than the RAM, things will be going in slow-motion.
The solution to that could be to either change the settings of the program you are running to decrease the memory usage (if the program has that functionality), or just get a computer with more memory.

Ok, the last command of the day (Yup, there will be more tomorrow!).
You might wonder how the heck am I supposed to be able to remember all these commands, options and flags? The simple answer is that you won't.
Not all of them at least.
You might remember ls, but was it -l or -a you should use to see hidden files? You might wish that there was a manual for these things.
Good news everyone, there is a manual! To get all the nitty-gritty details about ls, you use the **man** command.
Syntax:

```
$ man <command you want to look at>

$ man ls
```

```
LS(1)                            User Commands                            LS(1)

NAME
       ls - list directory contents

SYNOPSIS
       ls [OPTION]... [FILE]...

DESCRIPTION
       List  information  about  the  FILEs (the current directory by default).
       Sort entries alphabetically if none of -cftuvSUX nor  --sort  is  speci‐
       fied.

       Mandatory arguments to long options are mandatory for short options too.

       -a, --all
              do not ignore entries starting with .

       -A, --almost-all
              do not list implied . and ..

       --author
              with -l, print the author of each file

       -b, --escape
              print C-style escapes for nongraphic characters

       --block-size=SIZE
              scale  sizes by SIZE before printing them; e.g., '--block-size=M'
 Manual page ls(1) line 1 (press h for help or q to quit)
```

This will open a less window (remember, q to close it down, arrows to scroll) with the manual page about ls.
Here you will be able to read everything about ls.
You'll see which flag does what (-a is to show the hidden files, which in linux are files with a name starting with a dot .), which syntax the program has, etc.
If you are unsure about how to use a command, look it up in man.

The man pages can be a bit tricky to understand at first, but you get used to it with time.
If it is still unclear after looking in man, try searching for it on the internet.
You are bound to find someone with the exact same question as you, that has already asked on a forum, and gotten a good answer.
5 years ago.
The end.

If you still have time left on the lab and you finished early, check out the [advanced linux exercises](advanced-linux)!
