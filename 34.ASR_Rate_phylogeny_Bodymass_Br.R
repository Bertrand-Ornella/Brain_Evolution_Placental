########################  ASR under the variable rates model - Body mass  ########################

setwd("~/Desktop/Placental_April_2021")

library(phytools)
library(geoscale)
library(strap)

#open rates
Results1<-read.csv("Results_Body_Br_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#open calibrated tree
tree<-read.nexus("PEQ_tree.trees")
tree_scaled<-tree #need to be done because we will use tree to plot the ASR

##Rescale the branches of the timescaled phylogeny by mean scalar
scale<-Results$Median.Scalar/10000
tree_scaled$edge.length<-tree_scaled$edge.length*scale

#Estimate ancestral state reconstructions for this tree using ML 
Body.data1<-read.csv("Placental_data_final_PEQ.csv", header=T, row.names = 1)
Body.data <-Body.data1[6]
Body.data <-log10(Body.data)

Body.data<-as.matrix(Body.data)
Body<-as.numeric(Body.data)
Body<-setNames(Body, row.names(Body.data))
anc<-fastAnc(tree_scaled,Body,vars=TRUE,CI=TRUE)
anc2<-anc$ace

#save ancestral state reconstruction in excel file
num.anc2<-round(anc2,digits =2)
write.table(num.anc2, file="Body_Br_ancestral.csv")

#Plot these ancestral reconstructions in the regular time-scaled tree
plt<-fancyTree(tree,type="scattergram",X=as.matrix(Body),
               A=as.matrix(anc2),control=list(spin=FALSE),
               label="horizontal", plot=FALSE)

plt<-setMap(plt, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))
plt<-plt$contMaps[[1]]

#plot tree with node labels and time scale - !!! SAVE as A5 !!!
#Change "0.02" in ylim to change how tall the tree is
par(mfrow=c(1,1))
plot(plt, fsize=0.35,lwd=2.5,legend=FALSE,outline=FALSE,offset=1,
     ylim=c(-2-0.02*(Ntip(plt$tree)-1),Ntip(plt$tree)))

nodelabels(round(anc2,digits =2), cex=0.3, frame="none", adj = c(1.2, -0.75) )
tiplabels(round(Body[tree$tip.label],digits =2), cex=0.3, frame="none", adj= c(-0.4,0.3))

add.color.bar(leg=0.3*max(nodeHeights(tree)),plt$cols,title="Log10 (Body mass)",
              lims=plt$lims,digits=2,prompt=FALSE,lwd=5,fsize=0.6,subtitle="", x=5, y=25)

#Change "0.3" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

Tree.root.PEQ <-201.7734
axis(side=1, at=Tree.root.PEQ-15:Tree.root.PEQ, labels=15:Tree.root.PEQ,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#Now compare this ancestral reconstructions with reconstructions under a Brownian mode of evolution
anc_b<-fastAnc(tree,Body,vars=TRUE,CI=TRUE)
plb<-contMap(tree,Body, plot=FALSE)
plb<-setMap(plb, colors=c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1"))

plot(plb, fsize=0.4,lwd=2.5,outline=FALSE)
#nodelabels(round(anc_b$ace,digits =2), cex=0.5, frame="none", adj = c(1.2, -0.75) )

#To see differences in ancestral state reconstructions due the differences in the model assumptions 
anc2-anc_b$ace

#########################    Rate of evolution Body mass    ##############################
#open rates
Results1<-read.csv("Results_Body_Br_placentals.csv", header=T, row.names = 1)
Results<-Results1[-c(1),]

#Select data and log10 Median Scalar
Results$Median.Scalar<-log10(Results$Median.Scalar)
names(Results)[names(Results) == "Median.Scalar"] <- "Log.Median.Scalar"

Body.rate<-Results$Log.Median.Scalar

#open calibrated tree
tree<-read.nexus("PEQ_tree.trees")
par(mfrow=c(1,1))
plot(tree, cex=.3)

#Make plot
plotBranchbyTrait2<-function (tree, x, mode = c("edges", "tips", "nodes"), palette = "rainbow", 
                              legend = TRUE, xlims = NULL, ...) 
{
  mode <- mode[1]
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if (mode == "tips") {
    x <- c(x[tree$tip.label], fastAnc(tree, x))
    names(x)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  else if (mode == "nodes") {
    XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
    x <- rowMeans(XX)
  }
  if (hasArg(tol)) 
    tol <- list(...)$tol
  else tol <- 1e-06
  if (hasArg(prompt)) 
    prompt <- list(...)$prompt
  else prompt <- FALSE
  if (hasArg(type)) 
    type <- list(...)$type
  else type <- "phylogram"
  if (hasArg(show.tip.label)) 
    show.tip.label <- list(...)$show.tip.label
  else show.tip.label <- TRUE
  if (hasArg(show.node.label)) 
    show.node.label <- list(...)$show.node.label
  else show.node.label <- FALSE
  if (hasArg(edge.width)) 
    edge.width <- list(...)$edge.width
  else edge.width <- 4
  if (hasArg(edge.lty)) 
    edge.lty <- list(...)$edge.lty
  else edge.lty <- 1
  if (hasArg(font)) 
    font <- list(...)$font
  else font <- 3
  if (hasArg(cex)) 
    cex <- list(...)$cex
  else cex <- par("cex")
  if (length(cex) == 1) 
    cex <- rep(cex, 2)
  if (hasArg(adj)) 
    adj <- list(...)$adj
  else adj <- NULL
  if (hasArg(srt)) 
    srt <- list(...)$srt
  else srt <- 0
  if (hasArg(no.margin)) 
    no.margin <- list(...)$no.margin
  else no.margin <- TRUE
  if (hasArg(root.edge)) 
    root.edge <- list(...)$root.edge
  else root.edge <- FALSE
  if (hasArg(label.offset)) 
    label.offset <- list(...)$label.offset
  else label.offset <- 0.01 * max(nodeHeights(tree))
  if (hasArg(underscore)) 
    underscore <- list(...)$underscore
  else underscore <- FALSE
  if (hasArg(x.lim)) 
    x.lim <- list(...)$x.lim
  else x.lim <- NULL
  if (hasArg(y.lim)) 
    y.lim <- list(...)$y.lim
  else y.lim <- if (legend && !prompt && type %in% c("phylogram", 
                                                     "cladogram")) 
    c(1 - 0.06 * length(tree$tip.label), length(tree$tip.label))
  else NULL
  if (hasArg(direction)) 
    direction <- list(...)$direction
  else direction <- "rightwards"
  if (hasArg(lab4ut)) 
    lab4ut <- list(...)$lab4ut
  else lab4ut <- NULL
  if (hasArg(tip.color)) 
    tip.color <- list(...)$tip.color
  else tip.color <- "black"
  if (hasArg(plot)) 
    plot <- list(...)$plot
  else plot <- TRUE
  if (hasArg(rotate.tree)) 
    rotate.tree <- list(...)$rotate.tree
  else rotate.tree <- 0
  if (hasArg(open.angle)) 
    open.angle <- list(...)$open.angle
  else open.angle <- 0
  if (is.function(palette)) 
    cols <- palette(n = 1000)
  else {
    if (palette == "heat.colors") 
      cols <- heat.colors(n = 1000)
    if (palette == "gray") 
      cols <- gray(1000:1/1000)
    if (palette == "rainbow") 
      cols <- rainbow(1000, start = 0.7, end = 0)
  }
  if (is.null(xlims)) 
    xlims <- range(x) + c(-tol, tol)
  breaks <- 0:1000/1000 * (xlims[2] - xlims[1]) + xlims[1]
  whichColor <- function(p, cols, breaks) {
    i <- 1
    while (p >= breaks[i] && p > breaks[i + 1]) i <- i + 
        1
    cols[i]
  }
  colors <- sapply(x, whichColor, cols = cols, breaks = breaks)
  par(lend = 2)
  xx <- plot.phylo(tree, type = type, show.tip.label = show.tip.label, 
                   show.node.label = show.node.label, edge.color = colors, 
                   edge.width = edge.width, edge.lty = edge.lty, font = font, 
                   cex = cex[1], adj = adj, srt = srt, no.margin = no.margin, 
                   root.edge = root.edge, label.offset = label.offset, underscore = underscore, 
                   x.lim = x.lim, y.lim = y.lim, direction = direction, 
                   lab4ut = lab4ut, tip.color = tip.color, plot = plot, 
                   rotate.tree = rotate.tree, open.angle = open.angle, lend = 2, 
                   new = FALSE)
  if (legend == TRUE && is.logical(legend)) 
    legend <- round(0.3 * max(nodeHeights(tree)), 2)
  if (legend) {
    if (hasArg(title)) 
      title <- list(...)$title
    else title <- "trait value"
    if (hasArg(digits)) 
      digits <- list(...)$digits
    else digits <- 1
    if (prompt) 
      add.color.bar(legend, cols, title, xlims, digits, 
                    prompt = TRUE, fsize = cex[2])
    else add.color.bar(legend, cols, title, xlims, digits, 
                       prompt = FALSE, x = par()$usr[1] + 0.1 * (par()$usr[2] - 
                                                                   par()$usr[1]), y = par()$usr[3] + 0.4 * (par()$usr[4] - 
                                                                                                              par()$usr[3]), fsize = cex[2])
  }
  invisible(xx)
}


plotBranchbyTrait2(tree,Body.rate,mode="edges", title="Body mass evolutionary rate",
                   palette=colorRampPalette(c("black","mediumblue","skyblue","orange","goldenrod1","lightgoldenrod1")),
                   edge.width=2.5, cex=0.3, prompt=FALSE)

#Change "0.2" in mgp to change position of numbers on the scale bar (lower = higher)
#Change tck to change the length of the bars

Tree.root.PEQ <-201.7734
axis(side=1, at=Tree.root.PEQ-15:Tree.root.PEQ, labels=15:Tree.root.PEQ,cex.axis=0.5, pos=-2,
     mgp=c(2, 0.2, 0), tck=-0.02)

#label branches
edgelabels(round(Body.rate,digits =2), cex=0.3, frame="none", adj = c(0.9, -0.7) )


### END