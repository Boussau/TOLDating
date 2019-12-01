library(TreeSim)
#library(BioGeoBEARS)


n<-200
lambda <- 0.9 #0.3
mu <- 0.09  #0.02
frac <-1.0
age <- 4.0
numbsim<-1

i <- 1

# Each extant species is included in final tree with probability frac
# (the tree has n species AFTER sampling):

tre<-sim.bd.taxa.age(n,numbsim,lambda,mu,frac,age)

numTips <- length(tre[[1]]$tip.label)
tre[[1]]$node.label<-(numTips:(numTips+tre[[1]]$Nnode))

extant<-drop.extinct(tre[[1]])

extantLeaves <- setdiff( tre[[1]]$tip.label, is.extinct(tre[[1]]))
extinct <- drop.tip(tre[[1]], extantLeaves)

print("EXTANT:")
print(extant)

#outputting the trees
write.tree(tre[[1]], file = paste("TestData/completeTree_", as.character(i), ".dnd", sep=""))
if (extant$Nnode >=1) {
  write.tree(extant, file = paste("TestData/extantTree_", as.character(i), ".dnd", sep=""))
}


