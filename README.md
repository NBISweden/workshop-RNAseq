# workshop-RNAseq [![gh-actions-build-status](https://github.com/nbisweden/workshop-RNAseq/workflows/build/badge.svg)](https://github.com/nbisweden/workshop-RNAseq/actions?workflow=build)

This repo contains the material for NBIS workshop **Analysis of RNA-Seq data**. The rendered view of this repo is available [here](https://nbisweden.github.io/workshop-RNAseq/).

## Repo organisation

The 3-day course source material is located on the *master* branch (default). The 2-day course material is located on the *twoday* branch. The rendered material is located on the *gh-pages* branch. For most part, one only needs to update content in master. Changes pushed to the *master* or *twoday* branches are automatically rendered to the *gh-pages* branch.

The first build can take around 30-40 mins depending on the number of R packages. Subsequent builds take about 2 minutes since caching is enabled. Caches are removed after 7 days of last access. A push after that will require a full rebuild.

For more details about repo organisation, updating and modifying this repo, check out the [template repo](https://github.com/royfrancis/workshop-template-rmd-ga).

## Contributing

To add or update contents of this repo (for collaborators), first clone the repo.

```
git clone https://github.com/nbisweden/workshop-RNAseq.git
```

Make changes/updates as needed. Add the changed files. Commit it. Then push the repo back.

```
git add .
git commit -m "I did this and that"
git push origin
```

If you are not added as a collaborator, first fork this repo to your account, then clone it locally, make changes, commit, push to your repo, then submit a pull request to this repo.

## Rendering

The website is automatically rendered by [GitHub Actions](https://help.github.com/en/actions) whenever a change is pushed. **DO NOT** push any rendered material such as `slide_topic.html`, `lab_topic.html` or supporting directories `slide_topic_files`, `lab_topic_files` etc to GitHub.

For local rendering, you need R installed on your system along with dependency packages listed under *packages_cran_repo* and *packages_bioc_repo* in `_site.yml`. Steps below are run in R.

Run `rmarkdown::render_site()` in the project directory. This renders all .Rmd and .md files to generate the HTML files and all other necessary files (including the assets, images and data directories) and moves them all into a directory specified under `output_dir` in **`_site.yml`**. Open `output_dir/index.html` to start. Remove this directory after use. **DO NOT** commit and push this output directory to GitHub.

You can also run `rmarkdown::render("bla.Rmd")` on individual .Rmd/.md files. This is a time-saver as the whole website need not be rendered just to preview this one file.

---

**2020** NBIS • SciLifeLab
