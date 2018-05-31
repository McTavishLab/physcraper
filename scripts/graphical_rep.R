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
tre2<- drop.tip(sqTree, '1_9_schweinfurtii')
mtree<-make.simmap(tre2, traits, model = matrix(c(0, 0.00000001, 0.00000001, 0), 2))
plot(mtree,cols,type="fan",fsize=0.7,ftype="i")
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=0.8*par()$usr[3],fsize=0.8)