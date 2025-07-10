# Hirundo genomic divergence landscape analysis
# phylogeny.R - plot phylogenetic tree estimates for the genus

### Working directory and dependencies----
setwd('~/Google Drive/My Drive/projects/hirundo_divergence_landscape/analysis/phylogeny/')

library(ape)
library(phytools)

### Concatenated RAxML tree----

outgroup <- c('Petrochelidon_pyrrhonota')
# Sample taxa info
taxa <- read.table("./raxml/sample.info",header=F)

hirundo.rax <- read.tree("./raxml/hirundo.snps.focal.phy.raxml.support.tree")
hirundo.rax$tip.label <- taxa$V2
hirundo.rax <- root(hirundo.rax, outgroup, edgelabel = T)
hirundo.rax <- drop.tip(hirundo.rax,outgroup)
hirundo.rax$edge.length <- NULL
hirundo.rax <- ladderize(hirundo.rax)
plot(hirundo.rax,cex = 0.7, no.margin = T, label.offset = 0.15,edge.width=1)
nodelabels(hirundo.rax$node.label,frame= "none", cex = 0.5, adj = -0.25)

# Check to see order of input taxa - run this before the actual plotting above to format the 'sample.info' file
for (i in hirundo.rax$tip.label){
  print(i)
}

### Coalescent SVDquartets tree----

outgroup <- c('Petrochelidon_pyrrhonota')
# Sample taxa info
taxa <- read.table('./svdq/sample.info')

hirundo.svd <- read.nexus("./svdq/hirundo_genus.svdq.tre")
hirundo.svd$tip.label <- taxa$V2
hirundo.svd$edge.length <- NULL
hirundo.svd <- root(hirundo.svd, outgroup, edgelabel = T)
hirundo.svd <- drop.tip(hirundo.svd, outgroup)
boots.hirundo.svd <- hirundo.svd$node.label
hirundo.svd <- ladderize(hirundo.svd)
plot(hirundo.svd, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)
nodelabels(boots.hirundo.svd, frame= "none", cex = 0.5, adj = -0.25)

# Check to see order of input taxa - run this before the actual plotting above to format the 'sample.info' file
for (i in hirundo.svd$tip.label){
  print(i)
}

# Plot with only one H. rustica tip
rustica <- c('Hirundo_rustica_savignii','Hirundo_rustica_transitiva','Hirundo_rustica_tytleri','Hirundo_rustica_gutturalis','Hirundo_rustica_erythrogaster')
hirundo.svd2 <- drop.tip(hirundo.svd, rustica)
boots.hirundo.svd2 <- hirundo.svd2$node.label
hirundo.svd2$tip.label[1] <- 'Hirundo rustica'
plot(hirundo.svd2, cex = 0.6, no.margin = T, label.offset = 0.1, node.depth = 2)
nodelabels(boots.hirundo.svd2, frame= "none", cex = 0.5, adj = -0.25)

### Plot trees together
par(mfrow=c(1,2))
plot(hirundo.rax,cex = 0.6, no.margin = T, label.offset = 0.15,edge.width=1, node.depth = 2)
nodelabels(hirundo.rax$node.label,frame= "none", cex = 0.5, adj = -0.25)
plot(hirundo.svd, cex = 0.6, no.margin = T, label.offset = 0.15, node.depth = 2)
nodelabels(boots.hirundo.svd, frame= "none", cex = 0.5, adj = -0.25)
# output at 4 x 10
