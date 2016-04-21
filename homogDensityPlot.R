
# For graphing frequency data
setwd("/Users/christophermcallester/Documents/Personal Organization/Research/Thesis Code/TFCluster")
homogeneities = read.table("homogDensity.txt", sep = ",", header = T)

homDensity = density(homogeneities$Clustered, freq = F)
randDensity = density(homogeneities$Randomized, freq = F)

plot(randDensity,col=2,main="Density Plot of Homogeneity in Clusters and Random Reclusters",xlab="Homogeneity",ylab="Density")
points(homDensity,type='l')
legend("topright",c("Clustered","Randomized"), pch="",lty=1,col=c(1,2),bty="n")
