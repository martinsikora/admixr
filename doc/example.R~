
######################
## example analysis ##
######################

## --------------------------------------------------------------------------------
## load packages

library("admixr")
library("ggplot2")

## optional, for parallel computations
library("foreach")
library("doMC")


## --------------------------------------------------------------------------------
## load example data

basePath <- system.file(package="admixr")
extPath <- file.path(basePath, "extdata")
load(file.path(extPath, "ms.gt.RData"))
load(file.path(extPath, "ms.snpInfo.RData"))
load(file.path(extPath, "ms.sampleInfo.RData"))
load(file.path(extPath, "ms.popInfo.RData"))


## --------------------------------------------------------------------------------
## do PCA on the genotype matrix

res <- list()
res$full <- getPcaGT(gt) ## full PCA
res$apx <- getPcaGTFast(gt) ## faster approximate PCA

## plot
idx <- matrix(paste("PC", 1:10, sep = ""), nrow = 2)
idxC <- popInfo$color ## colors for populations
names(idxC) <- popInfo$popId

for(k in names(res)){
    d <- as.data.frame(res[[k]]$summary$pca)
    d$sampleId <- rownames(d)
    d$popId <- sampleInfo$popId[match(d$sampleId, sampleInfo$sampleId)]
    pdf(paste(k, "pca.plot.pdf", sep = "."), width = 4, height = 4)
    for(i in 1:ncol(idx)){
        p <- ggplot(d, aes_string(x = idx[1, i], y = idx[2, i], colour = "popId", fill = "popId"))
        print(p + geom_text(aes(label = gsub("pop", "", popId)), size = 2, alpha = 0.5) + scale_colour_manual(name = "Population", values = idxC) + scale_fill_manual(name = "Population", values = idxC) + labs(x = idx[1, i], y = idx[2, i]) + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title = element_text(size = 10), legend.position = "none"))
    }
    dev.off()
}


## --------------------------------------------------------------------------------
## get allele frequencies

idxP <- sampleInfo$popId[match(colnames(gt), sampleInfo$sampleId)] ## index vector for population labels of individuals
counts <- getAlleleCountsGT(gt, idxP) ## allele counts for both alleles in each population
freq <- counts[,,2] / (counts[,,1] + counts[,,2]) ## frequency of allele 2 in each population


## --------------------------------------------------------------------------------
## calculate FST for all pairs of populations

idxP <- sampleInfo$popId[match(colnames(gt), sampleInfo$sampleId)] ## index vector for population labels of individuals
idxM <- t(combn(popInfo$popId, 2)) ## index matrix of all pairwise comparisons
fst <- as.data.frame(idxM, stringsAsFactors = FALSE)
colnames(fst) <- c("p1", "p2")
fst$fst <- NA
fst1 <- fst ## second data frame for calculation from allele counts

for(i in 1:nrow(idxM)){
    cat(i, "\r")
    idx <- idxP %in% idxM[i,]
    fst$fst[i] <- getFstGT(gt[, idx], idxP[idx], region = TRUE) ## Weir & Cockerham 1984, for diploid genotypes
}
for(i in 1:nrow(idxM)){
    cat(i, "\r")
    fst1$fst[i] <- getFstAlleleCounts(counts[, idxM[i,],], region = TRUE) ## Weir & Hill 2002, from allele frequencies assuming no inbreeding
}

## compare the estimators
pdf("fst.comparison.plot.pdf", width = 5, height = 5)
plot(fst$fst, fst1$fst, xlab = "WC1984", ylab = "WH2002", asp = 1, type = "n")
abline(0, 1, col = "grey")
points(fst$fst, fst1$fst, pch = 21, col = "#00000055", bg = "#00000055")
dev.off()

## convert to matrix & make nj tree from WC1984 estimator results
fst1 <- fst[, c(2:1, 3)]
colnames(fst1) <- colnames(fst)
fst1 <- rbind(fst, fst1)
fst1$p1 <- factor(fst1$p1, levels = popInfo$popId)
fst1$p2 <- factor(fst1$p2, levels = popInfo$popId)
d <- reshape2::acast(fst1, p1 ~ p2, value.var = "fst")
d1 <- ape::bionj(d) ## bionj algorithm from ape package

## plot
pdf("fst.nj.plot.pdf", width = 8, height = 7)
plot(ape::root(d1, outgroup = "pop1"), font = 1, tip.color = idxC[d1$tip.label], cex = 1)
plot(d1, type = "unrooted", lab4ut = "axial", font = 1, tip.color = idxC[d1$tip.label], cex = 1)
dev.off()

pal <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
pdf("fst.matrix.plot.pdf", width = 6, height = 5)
p <- ggplot(fst1, aes(x = p1, y = p2, fill = fst))
print(p + geom_tile() + geom_text(aes(label = formatC(fst, digits = 2, format = "f")), size = 2) + scale_fill_gradientn(name = expression(F[ST]), colours = pal, limits = c(0, max(fst$fst))) + coord_equal() + xlab("Population 1") + ylab("Population 2") + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major = element_blank(), strip.background = element_rect(fill = "white", colour = NA), strip.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = idxC), axis.text.y = element_text(color = idxC), axis.ticks.length = grid::unit(0, "cm")))
dev.off()


## --------------------------------------------------------------------------------
## outgroup f3 statistic

blockIdx <- getChromosomeBlocks(snpInfo$chromosome, snpInfo$position, 5e6)
idxO <- popInfo$popId[1]
idxX <- c("pop10", "pop20")
f3 <- list()
registerDoMC(4) ## run 4 threads in parallel for f3 computation
for(x in idxX){
    idxP1 <- popInfo$popId[!(popInfo$popId %in% c(idxO, x))]
    f3[[x]] <- foreach(p1 = idxP1, .combine = "rbind") %dopar% {
        idxM <- cbind(idxO, x, p1)
        r <- doF3Test(counts, idxM, blockIdx)
        r
    }
}
f3 <- do.call("rbind", f3)

## plot
pdf("f3.outgroup.plot.pdf", width = 4, height = 3.5)
for(x in idxX){
    d <- f3[f3$p2 == x,]
    print(plotFStat(d$f3, d$se, d$p3, d$p3, idxC, intercept = max(d$f3), orderValue = TRUE, showLegend = FALSE) + ylab(expression(f[3])) + ggtitle(x)) ## add ylab and title to plot object returned from function
}
dev.off()

## plot pop 10 results grouped by outgroup / derived population
d <- f3[f3$p2 == "pop10",]
groups <- ifelse(as.integer(gsub("pop", "", d$p3)) < 10, "outgroup", "derived")
idxC1 <- c(outgroup = "blue", derived = "red")
pdf("f3.outgroup.grouped.plot.pdf", width = 4, height = 3.5)
print(plotFStat(d$f3, d$se, d$p3, groups, idxC1, intercept = max(d$f3), grouped = TRUE, showLegend = FALSE) + ylab(expression(f[3])) + ggtitle("pop10")) ## add ylab and title to plot object returned from function
dev.off()

## plot pairwise f3 for pop10 / pop20
d1 <- f3[f3$p2 == "pop10" & f3$p3 != "pop20",]
d2 <- f3[f3$p2 == "pop20" & f3$p3 != "pop10",]

pdf("f3.outgroup.pairs.plot.pdf", width = 4, height = 4)
print(plotFPairs(d1$f3, d2$f3, d1$se, d2$se, d1$p3, idxC, z = 1, shape = 21, alpha = 0.9, showLegend = FALSE) + xlab("pop10") + ylab("pop20")) ## add axis labels to plot object returned from function
print(plotFPairs(d1$f3, d2$f3, d1$se, d2$se, d1$p3, label = d1$p3, idxC, z = 1, size = 0, alpha = 0.9, showLegend = FALSE) + xlab("pop10") + ylab("pop20")) ## add axis labels to plot object returned from function
dev.off()


## --------------------------------------------------------------------------------
## D statistic

idxM <- cbind("pop1", "pop10", popInfo$popId[c(2:4, 6:9, 11:20)], "pop5")
registerDoMC(4)
d <- foreach(i = 1:nrow(idxM), .combine = "rbind") %dopar% {
    r <- doDTest(freq, idxM[i, , drop = FALSE], blockIdx)
    r
}

## plot
pdf("D.plot.pdf", width = 3, height = 4)
print(plotFStat(d$D, d$se, d$p3, d$p3, idxC, horizontal = FALSE, showLegend = FALSE) + xlab("D") + ggtitle("D(pop1,pop10)(X,pop5)")) ## add xlab and title to plot object returned from function
dev.off()


## --------------------------------------------------------------------------------
## ADMIXTURE projections of pop10 and pop20

## read full panel results
qMatrix <- matrix(scan(file.path(extPath, "ms.adm.5.Q")), ncol = 5, byrow = TRUE)
fam <- read.table(file.path(extPath, "ms.adm.fam"), stringsAsFactors = FALSE)
rownames(qMatrix) <- fam$V2

## plot
d <- reshape2::melt(qMatrix)
colnames(d) <- c("sampleId", "k", "value")
d$sampleId <- factor(d$sampleId, levels = unique(d$sampleId))
d$popId <- sampleInfo$popId[match(d$sampleId, sampleInfo$sampleId)]
d$popId <- factor(d$popId, levels = unique(d$popId))
idxCLabel <- popInfo$color[match(levels(d$popId), popInfo$popId)]
idxCCluster <- c("violetred3", "wheat3", "steelblue3", "paleturquoise3", "gold2")
pdf("admixture.full.plot.pdf", width = 8, height = 2.5)
plotAdmixture(sampleId = d$sampleId, popId = d$popId, k = d$k, value = d$value, colors = idxCCluster, labColors = idxCLabel, alpha = 1, width = 1, showLegend = TRUE, rot = 90)
dev.off()


## read ref panel without pop10, pop20 results
qMatrix <- matrix(scan(file.path(extPath, "ms.adm.ref.5.Q")), ncol = 5, byrow = TRUE)
fam <- read.table(file.path(extPath, "ms.adm.ref.fam"), stringsAsFactors = FALSE)
rownames(qMatrix) <- fam$V2

## prepare data for individuals to be projected
idxS <- sampleInfo$sampleId[sampleInfo$popId %in% c("pop10", "pop20")]
bim <- read.table(file.path(extPath, "ms.adm.ref.bim"), stringsAsFactors = FALSE)
idx <- match(bim$V2, snpInfo$posId) ## find indices of SNPs used in ref panel
gt1 <- gt[idx, idxS] ## subset of genotype matrix for individuals to project and ref panel SNPs

idx <- bim$V5 == 1 ## find SNPs where alleles are switched in recoded subset plink file of ref panel; genotype matrix has copies of minor allele, which is coded as allele "2" in the plink file
gt1[idx,] <- 2 - gt1[idx,]
gt1 <- 2 - gt1 ## final step is to reverse all snps for estimation, to match coding of P matrix from admixture (has major allele instead of minor allele coding)

## do projection
pMatrix <- matrix(scan(file.path(extPath, "ms.adm.ref.5.P")), ncol = 5, byrow = TRUE)
registerDoMC(4) ## run 4 threads in parallel
qMProj <- foreach(s = idxS, .combine = "rbind") %dopar% {
    cat(s, "\r")
    doAdmixtureProjection(gt1[, s], pMatrix)
}
rownames(qMProj) <- idxS
qMatrix <- rbind(qMatrix, qMProj)
qMatrix <- qMatrix[sampleInfo$sampleId,]

## plot
d <- reshape2::melt(qMatrix)
colnames(d) <- c("sampleId", "k", "value")
d$sampleId <- factor(d$sampleId, levels = unique(d$sampleId))
d$popId <- sampleInfo$popId[match(d$sampleId, sampleInfo$sampleId)]
d$popId <- factor(d$popId, levels = unique(d$popId))
idxCLabel <- popInfo$color[match(levels(d$popId), popInfo$popId)]
idxCCluster <- c("violetred3", "steelblue3", "paleturquoise3", "wheat3", "gold2")
pdf("admixture.projected.plot.pdf", width = 8, height = 2.5)
plotAdmixture(sampleId = d$sampleId, popId = d$popId, k = d$k, value = d$value, colors = idxCCluster, labColors = idxCLabel, alpha = 1, width = 1, showLegend = TRUE, rot = 90)
dev.off()
