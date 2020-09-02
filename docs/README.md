# Compiling docs

Build vignettes:
```
rmarkdown::render('vignettes/spatially_aware_clustering.Rmd', clean=FALSE, output_format='all')
```

Move files from vignettes directory to docs:
```
# delete old files
rm -rf spatially_aware_clustering_files/
mv ../vignettes/spatially_aware_clustering_files/ .
# overwrite old files
mv ../vignettes/spatially_aware_clustering.knit.md  spatially_aware_clustering.md
mv ../vignettes/spatially_aware_clustering.pdf .
```