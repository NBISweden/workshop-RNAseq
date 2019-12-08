#!/bin/bash

## this script renders Rmd/md files into a new directory
## pulls down gh-pages branch and copies the new directory into it
## then pushes it to the course repo's gh-pages branch

set -e

# get repo url, create new url with token
url_git=$(git config --get remote.origin.url)
echo "Repo: $url_git"
url_git_short=$(echo $url_git | sed 's/https:\/\///')
url_repo="https://${GITHUB_TOKEN}@${url_git_short}"

# get current directory and destination directory
path_root=$(pwd)
path_export=$(grep -E '^output_dir' _site.yml | sed -E 's/^[a-zA-Z0-9_-]+[:][[:space:]]//')
echo "Current directory is ${path_root}"

# render html content usin R if _site.yml file is available else just copy
if [ -f _site.yml ]; then
  echo "Rendering content ..."
  R -q -e "rmarkdown::render_site()"
else
  echo "Copying contents ..."
  mkdir $path_export
  find . ! -name $path_export -exec cp -r -t $path_export {} +
fi

# clone gh-pages branch alone from repo
git clone --single-branch --branch gh-pages $url_repo tmprepo
cd tmprepo
pwd

# remove destination directory if it exists and create new
if [ -d $path_export ]; then
  echo "Directory ${path_export} already exists. Removing the directory."
  git rm -r $path_export
  git commit -m "Old directory ${path_export} deleted."
fi

cd ${path_root}
cp -r $path_export tmprepo/
cd tmprepo

echo "Creating index.md ..."
# create index.md
str_user_repo=$(echo $url_git | sed 's/https:\/\/github.com\///' | sed 's/.git//')
str_user=$(echo "$str_user_repo" | sed 's/\/.*//')
str_repo=$(echo "$str_user_repo" | sed 's/.*\///')

printf "All current and previous workshop materials are listed below. Numbers denote year and month.\n\n" > index.md

dirs=$(ls -d */ | sed 's/\///' | tac)

for i in $dirs
do
 printf "<a href='https://${str_user}.github.io/${str_repo}/${i}/'>${i}</a>\n\n" >> index.md
done

echo "Deploying to gh-pages ..."
git add .
tm=$(date +%Y%m%d-%H%M%S)
git commit -m "Updated contents at ${tm}"
git push origin

cd ${path_root}
rm -rf tmprepo

echo "Completed successfully."
