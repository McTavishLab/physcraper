library(ape)
library(phytools)
sqData<-read.csv("output.csv", stringsAsFactors = FALSE)
sqTree<-read.newick("physcraper.tre")
plotTree(sqTree,type="fan",lwd=1,fsize=0.5)
library(geiger)
traits <- sqData[,2]
species<-sub("'","",sqData[,1])
species<-sub("'","",species)
names(traits)<-species
fitER<-fitDiscrete(sqTree, traits,model="ER")

mtree<-make.simmap(sqTree,traits,model="ER")
mtree<-make.simmap(sqTree,traits,model = matrix(c(0, 0.00000001, 0.00000001, 0), 2))
plot(fitER)
sqTree.label
sqTree
