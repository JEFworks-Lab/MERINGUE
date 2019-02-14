# Compiling docs

Build vignettes:
```
rmarkdown::render('vignettes/mOB_analysis.Rmd', clean=FALSE, output_format='all')
```

Move files from vignettes directory to docs:
```
# delete old files
rm -r mOB_analysis_files/
mv ../vignettes/mOB_analysis_files/ .
# overwrite old files
mv ../vignettes/mOB_analysis.md  .
mv ../vignettes/mOB_analysis.pdf .
```