---
layout: default
title: Upload R package to github
parent: Resources
nav_order: 3
permalink: /docs/Resources/Rpackage
---

##### First make an R package:

```
library(Rcpp)
library(RcppArmadillo)
RcppArmadillo.package.skeleton("mypackage", example_code = FALSE)
setwd("mypackage")
# delete Read-and-delete-me
# Modify Description file
# Copy the R codes to folder R
# Copy the Cpp codes to folder src
usethis::use_gpl3_license("Lulu Shang") # Add license 
compileAttributes(verbose=TRUE) # Find and register Rcpp functions 
devtools::load_all() # load all functions
# pkgbuild::compile_dll()
devtools::document() # create roxygen2 document, don't make .md by hand, let roxygen2 do it
devtools::check() # check whether the package is OK
devtools::build() # build a package
```
##### Push R package to github:

Go to the R package dictionary, init the git
```
cd biostat815hw1/
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git init
Initialized empty Git repository in /home/lulushang/Projects/course/hyun815/hw1/biostat815hw1/.git/

lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git add .
```

```
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git config --global user.email "youremail@xxx.edu"
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git config --global user.name "yourname"
```

Now we can write the first commit
```
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git commit -m "Initial commit"
[master (root-commit) 74f4fc1] Initial commit
 17 files changed, 866 insertions(+)
 create mode 100644 .Rbuildignore
 create mode 100644 DESCRIPTION
 create mode 100644 LICENSE.md
 create mode 100644 NAMESPACE
 create mode 100644 R/RcppExports.R
 create mode 100644 R/logisticLLKr.r
 create mode 100644 R/logisticNelderMead.r
 create mode 100644 man/biostat815hw1-package.Rd
 create mode 100644 man/logisticLLKr.Rd
 create mode 100644 man/logisticNelderMead.Rd
 create mode 100644 src/Makevars
 create mode 100644 src/Makevars.win
 create mode 100644 src/RcppExports.cpp
 create mode 100644 src/RcppExports.o
 create mode 100755 src/biostat815hw1.so
 create mode 100644 src/logisticLLKc.cpp
 create mode 100644 src/logisticLLKc.o
```

Remote add origin
```
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$  git remote add origin https://github.com/shangll123/biostat815hw1
```

Try to pull 
```
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git pull --allow-unrelated-histories
Username for 'https://github.com': shangll123
Password for 'https://shangll123@github.com': 
There is no tracking information for the current branch.
Please specify which branch you want to merge with.
See git-pull(1) for details.

    git pull <remote> <branch>

If you wish to set tracking information for this branch you can do so with:

    git branch --set-upstream-to=origin/<branch> master
```

Finish
```
lulushang@master:~/Projects/course/hyun815/hw1/biostat815hw1$ git push --set-upstream origin master
Counting objects: 21, done.
Delta compression using up to 24 threads.
Compressing objects: 100% (20/20), done.
Writing objects: 100% (21/21), 879.01 KiB | 3.50 MiB/s, done.
Total 21 (delta 1), reused 0 (delta 0)
remote: Resolving deltas: 100% (1/1), done.
To https://github.com/shangll123/biostat815hw1
 * [new branch]      master -> master
Branch 'master' set up to track remote branch 'master' from 'origin'.


```

Useful link: https://stackoverflow.com/questions/10298291/cannot-push-to-github-keeps-saying-need-merge


### Additional info:
After you have created an R package, you can directly decompress the zipped file and upload the folder to your github repository. This might be the easiest way as far as I know. 

