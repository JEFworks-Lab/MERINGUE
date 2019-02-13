# test different types of nearest neighbors

source('../R/helper.R')

ct <- as.factor(c(1,1,1,2,2,2))
pos <- data.frame(x = c(1,1,2,3,3,2),
                  y = c(1,2,1,3,4,4))
names(ct) <- rownames(pos) <- paste0('cell', 1:6)
ctA <- names(ct)[ct==1]
ctB <- names(ct)[ct==2]

plot(pos, col=rainbow(2)[ct], pch=16, cex=2)

n1 <- getAdj(pos, k=2)
plotNetwork(pos, n1, col=rainbow(2)[ct], main='Nearest Neighbors (k=2)', cex=2)

n2 <- getMnn(ctA, ctB, pos, k=2)
plotNetwork(pos, n2, col=rainbow(2)[ct], main='Mutual Nearest Neighbors (k=2)', cex=2)

n3 <- getBnn(ctA, ctB, pos, k=2)
plotNetwork(pos, n3, col=rainbow(2)[ct], main='Background Nearest Neighbors (k=2)', cex=2)

n4 <- getBnn(ctB, ctA, pos, k=2)
plotNetwork(pos, n4, col=rainbow(2)[ct], main='Background Nearest Neighbors (k=2)', cex=2)
