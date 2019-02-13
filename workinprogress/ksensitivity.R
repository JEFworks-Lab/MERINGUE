# Test sensitivity to k

cd <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/data/Rep11_MOB_0.csv', row.names=1)
pos <- read.csv('/Users/jefworks/Desktop/SpatialDE/Analysis/MouseOB/MOB_sample_info.csv', row.names=1)
cd <- cd[rownames(pos),]

results.all <- lapply(3:10, function(k) {
  print(k)
  # plot
  adj <- getAdj(pos[,1:2], k=k)
  plotNetwork(pos[,1:2], adj)

  # calculate Moran's I
  results <- do.call(rbind, parallel::mclapply(seq_len(nrow(mat)), function(i) {
    value <- mat[i,]
    ape::Moran.I(value, adj)
  }, mc.cores=parallel::detectCores()-1))
  rownames(results) <- rownames(mat)
  results <- as.data.frame(results)

  results$p.adj <- stats::p.adjust(results$p.value)

  return(results)
})

plot(-log10(unlist(results.all[[1]]$p.value)), -log10(unlist(results.all[[5]]$p.value)))
