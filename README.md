# workshop-RNAseq [![gh-actions-build-status](https://github.com/nbisweden/workshop-RNAseq/workflows/build/badge.svg)](https://github.com/nbisweden/workshop-RNAseq/actions?workflow=build)

This repo contains the material for NBIS workshop **Analysis of RNA-Seq data**. The rendered view of this repo is available [here](https://nbisweden.github.io/workshop-RNAseq/).

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

:exclamation: When updating repo for a new course, change `output_dir: XXXX` in `_site.yml` as the first thing, so that old rendered files are not overwritten.

:exclamation: Do not push any rendered .html files or intermediates.

### Local build/preview using Docker

You can preview changes and build the whole website locally without a local installation of R or dependency packages by using the pre-built Docker image.

:exclamation: **Note:** Large image size: 4.6GB.

Clone the repo if not already done. Make sure you are standing in the repo directory.

To build the complete site,

```
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${PWD}:/rmd ghcr.io/nbisweden/workshop-rnaseq:latest
```

To build a single file (for example `index.Rmd`),

```
docker run --rm -u $(id -u ${USER}):$(id -g ${USER}) -v ${PWD}:/rmd ghcr.io/nbisweden/workshop-rnaseq:latest Rscript -e 'rmarkdown::render("index.Rmd")'
```

See **Dockerfile** to build the image.

:exclamation: Output files are for local preview only. Do not push any rendered .html files or intermediates.

## Repo organisation

The 3-day course source material is located on the *master* branch (default). The 2-day course material is located on the *twoday* branch. The rendered material is located on the *gh-pages* branch. For most part, one only needs to update content in master. Changes pushed to the *master* or *twoday* branches are automatically rendered to the *gh-pages* branch.

:exclamation: Every push rebuilds the whole website using a docker image. Build takes about 6 mins.


---

**2024** • NBIS • SciLifeLab
