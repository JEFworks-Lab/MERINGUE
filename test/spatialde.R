library(reticulate)
use_python("/usr/bin/python")

## Make sure the proper python path is loaded where SpatialDE is installed
# PYTHONPATH="/usr/local/lib/python2.7/site-packages/":"${PYTHONPATH}"
# export PYTHONPATH

SpatialDE <- import("SpatialDE")

umap_res <- umap$UMAP()$fit_transform(X)
results = SpatialDE.run(pos, mat)
