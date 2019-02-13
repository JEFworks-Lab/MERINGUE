## test bplapply
library(BiocParallel)

fun <- function(v) {
    setTxtProgressBar(pb, v)
    message("working") ## 10 tasks
    Sys.sleep(1)
}

start = Sys.time()
ncores = 10
pb <- txtProgressBar(min=0, max=100, char = ".") # use progress bar
BiocParallel::bplapply(1:100, fun, BPPARAM = BiocParallel::MulticoreParam(workers=ncores, tasks=1, progressbar=TRUE))
close(pb)
end = Sys.time()
end - start

start = Sys.time()
ncores = 1
pb <- txtProgressBar(min=0, max=10, char = ".") # use progress bar
BiocParallel::bplapply(1:10, fun, BPPARAM = BiocParallel::MulticoreParam(workers=ncores, progressbar=TRUE))
close(pb)
end = Sys.time()
end - start


myFun = function(i) { Sys.sleep(runif(1) / 4); message(i) }
xx = bplapply(1:10, myFun, BPPARAM=MulticoreParam(4, tasks=4))
xx = bplapply(1:10, myFun, BPPARAM=MulticoreParam(4, tasks=1))
xx = bplapply(1:100, myFun, BPPARAM=MulticoreParam(4, tasks=100, progressbar=TRUE))
xx = bplapply(1:100, myFun, BPPARAM=MulticoreParam(4, tasks=1, progressbar=TRUE))
xx = bplapply(1:100, myFun, BPPARAM=MulticoreParam(4, tasks=10, progressbar=TRUE))
