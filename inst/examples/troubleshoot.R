                                                                     
                                                                     
                                                                     
                                             
# load libraries
library(ape)
library(ouch)

# load tree
tree<-read.nexus("ttree.nex")

# load data
data<-read.table("data_ou.txt")

# prepare the tree and the data

# get tree into ouch format
outree<-ape2ouch(tree)

# convert the tree into a dataframe
oudf<-as(outree,"data.frame")

# add the taxon names to a new column in the continuous dataset called labels
data$labels<-row.names(data)

# merge the dataset with the oudataframe using the labels columns in both objects
oudata<-merge(oudf,data,by="labels",all=T)

# set the row names in this dataframe to the node labels
row.names(oudata)<-oudata$nodes

# remake the outree from the ground up from the ou dataframe we created
newoutree<-ouchtree(nodes=oudata$nodes, ancestors=oudata$ancestors, times=oudata$times, labels=oudata$labels)

# check tree
plot(newoutree, node.names=T)

# assign each of the nodes and tips to an optima

# ou1
oudata$ou1<-as.factor("single")

# X1
X1<-oudata$X1
X1[is.na(X1)==T]<-c("S","G","G","G","G","G","G","S","G","S","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","S","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","G","S","S","G","G","G","G","G","S","S","S","S","S","S")
oudata$X1<-as.factor(X1)

# X2
X2<-oudata$X2
X2[is.na(X2)==T]<-c("A","O","O","O","O","O","O","P","O","P","O","O","O","O","O","O","O","O","O","O","O","O","O","O","O","O","P","O","O","O","O","O","O","O","O","O","O","O","O","O","A","O","A","O","A","O","O","A","A","O","O","O","O","O","A","A","A","A","A","A")
oudata$X2<-as.factor(X2)

# X3
X3<-oudata$X3
X3[is.na(X3)==T]<-c(1,3,3,3,3,3,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,2,2,2,2,3,3,1,1,1,1,2,2,1,1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1,1)
oudata$X3<-as.factor(X3)

# X4
X4<-oudata$X4
X4[is.na(X4)==T]<-c(1,3,3,3,3,3,3,4,3,4,3,3,3,3,3,3,3,3,3,3,3,3,3,1,1,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,2,2,2,2,2,1,1,1,1,1,1)
oudata$X4<-as.factor(X4)

# load packages
library(pmc)

# load parallelization package
library(snowfall)
sfInit(parallel=T,cpu=4)
sfLibrary(pmc)

# load auxiliary packages
library(reshape2)
library(reshape)
library(grid)

# compare BM vs Pagel
evr.0_lambda<-pmc(newoutree,oudata["evr"],modelA=c("brown"),modelB=c("lambda"),optionsA=list(),optionsB=list(),nboot=5000)
a.0_lambda<-evr.0_lambda[["A"]]
b.0_lambda<-evr.0_lambda[["B"]]
evr.0v1.pd<-evr.0v1$par_dists
write.table(evr.0_lambda.pd,"0_lambda-par_dists.txt",row.names=T,quote=F)
evr.0_lambda.ic<-cast(evr.0_lambda$par_dist, comparison ~ parameter,function(x) quantile(x, c(0.05, 0.95)), value = c("lower","upper"))
write.table(evr.0_lambda.ic,"0_lambda-ic.txt",row.names=T,quote=F)
