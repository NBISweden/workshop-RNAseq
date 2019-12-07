# course_rnaseq

This repo contains the course material for NBIS course **Workshop on RNA-Seq**. The rendered view of this repo is available [here](https://nbisweden.github.io/course_rnaseq/).

## Contributing

To add or update contents of this repo (for collaborators), first clone the repo.

```
git clone https://github.com/nbisweden/course_rnaseq.git
```

Make changes/updates as needed. Add the changed files. Commit it. Then push the repo back.

```
git add *
git commit -m "I did this and that"
git status
git push origin
```

If you are not added as a collaborator, first fork this repo to your account, then clone it locally, make changes, commit, push to your repo, then submit a pull request to this repo.

## Descriptions

These are descriptions of the files and a guide to updating this repo for course maintainers.

**_site.yml**  
All website configuration options are here. Nothing to be changed here.

**README.md**  
You are reading this file now. Nothing to do here.

**index.Rmd**  
This file generates the home page. **Make sure that the date and location is correct**. Verify links.

**schedule.Rmd**  
This file generates the schedule page. **The start time of the course is set to `09:00:00`. Change if needed**. It is assumed that the course starts at the same time everyday.  

**schedule.csv**  
This table holds the schedule information. Open/edit in a spreadsheet or text editor. Columns are delimited by `;`. Do not change the number of columns, position of columns or column names.

***date***: Full date for each day in format dd/mm/yyyy. Missing/empty cells are filled down automatically.  
***room***: Room number for the course. Missing/empty cells are filled down automatically.  
***dur***: Duration for the topic in minutes.  
***topic***: Topic name (Keep it short).  
***person***: Name of the person covering the topic.  
***link_presentation***: (Optional) Link to the presentation. Local links can be just `presentation_topic.html`. Use this labelling convention.  
***link_lab***: (Optional) Link to the lab material. Local links can be just `lab_topic.html`. This is the labelling convention used.  
***link_room***: (Optional) Link to the room location. Can be a google map link, mazemap link etc.  

**lab.Rmd**  
This page brings together all the exercises. Verify links.

**precourse.Rmd**  
This page holds the steps needed to be completed before the course. **Change project ID for each course**.

**info.Rmd**  
This page contains practical information related to the course. **Set location and update info if needed**.

[**presentation_topic.Rmd**]  
RMarkdown presentation files for various topics. Replace 'topic' with a short name of the topic.  

[**lab_topic.Rmd**]  
RMarkdown lab files for various topics. Replace 'topic' with a short name of the topic.  

## Dependencies

The `assets` directory contains css styles, headers, footers, logos etc. If you are using images in your .Rmd file, place them in the directory `images` and refer to them using relative path like `![](./images/image.jpg)`. Images generated in R during rendering of the .Rmd file is automatically handled. If you have data (tsv, csv, txt text files, .Rds files), place them inside the directory `data` and read them using relative path `x <- read.delim("./data/table.txt")`. Do not use paths that link outside of the project environment.

## Rendering

The website is rendered by running `rmarkdown::render_site()` in the project directory. This generates the HTML files and all other necessary files (including the assets, images and data directories) and moves them into a directory named **docs**. Open **docs/index.html** to start. The output directory is set to **docs** because this GitHub repo uses the **docs** directory as the Github pages (rendered content) source.

For testing purposes, you can run `rmarkdown::render("bla.Rmd")` on individual Rmd files. This is a time-saver as the whole website need not be rendered just to preview this one file.

---

**2019** NBIS | SciLifeLab
# rnaseq_201911
