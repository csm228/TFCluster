
# For graphing frequency data
setwd("/Users/christophermcallester/Documents/Personal Organization/Research/Thesis Code/TFCluster")
homogeneities = read.table("percentHomogeneities.txt", fill = TRUE, stringsAsFactors = FALSE)

nams <- homogeneities[, 1]
homogeneities <- homogeneities[, -1]
listData <- split(homogeneities, seq_len(nrow(homogeneities)))
listData <- lapply(listData, function(x) x[x != ""])
names(listData) <- nams

print(listData)

plot(listData$Dataset10,ylim=c(0,1.0),col='red',type='l',main="Percent Homogeneity of Clusters in Ten Datasets",xlab="Cluster",ylab="Homogeneity")
points(listData$Dataset9,type='l',col='orange')
points(listData$Dataset8,type='l',col='yellow')
points(listData$Dataset7,type='l',col='green')
points(listData$Dataset6,type='l',col='blue')
points(listData$Dataset5,type='l',col='lightblue')
points(listData$Dataset4,type='l',col='darkblue')
points(listData$Dataset3,type='l',col='violet')
points(listData$Dataset2,type='l',col='purple')
points(listData$Dataset1,type='l',col='lightgreen')
legend("topright",c("1","2","3","4","5","6","7","8","9","10","Null"), pch="",lty=1,col=c('lightgreen','purple','violet','darkblue','lightblue','blue','green','yellow','orange','red','black'),bty="n")

