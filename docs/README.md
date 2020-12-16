# Compiling docs

Build vignettes:
```
rmarkdown::render('vignettes/simulation.Rmd', clean=FALSE, output_format='md_document')
```

Move files from vignettes directory to docs:
```
# delete old files
rm -rf simulation_files/
mv ../vignettes/simulation_files/ .
# overwrite old files
#mv ../vignettes/simulation.knit.md  simulation.md
mv ../vignettes/simulation.md  .
#mv ../vignettes/simulation.pdf .
```