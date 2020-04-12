
holder="/tmp/tmp.1rvR7T5LoE"

for(param in c("age","IMPRES","IFNg","Bindea","ESTIMATE","HRDscore","SNVs","indels","dels","ins"))
{
  d=read.table(paste0(holder,"_",param,"_summary"))

  box_source=boxplot(d,plot=F)
  #change whiskers to 5 and 95th percentile
  for(i in c(1,4))
  {
    box_source$stats[1,i]<-quantile(d[,i],probs=c(0.05,0.95),na.rm=TRUE)[1]
    box_source$stats[5,i]<-quantile(d[,i],probs=c(0.05,0.95),na.rm=TRUE)[2]
  }

  pdf(paste0("~/R_output/",param,"_boxplot.pdf"),useDingbats=F)
  bxp(box_source,frame.plot=F,ylab=param,main=param,boxfill=c((rainbow(4))),
  names=colnames(d))

  for(i in c(1,2,3,4))
  {
    line=paste0("median=",box_source$stats[3,i],
    "\nn=",length(which(!is.na(d[,i]))))
    text(i,box_source$stats[3,i],line)
  }

  dev.off()

}

#post-test special cases
#SNVs, Clusters 2 vs 4
param="SNVs"
d=read.table(paste0(holder,"_",param,"_summary"))
wilcox.test(d$cluster2,d$cluster4)
W = 2609, p-value = 8.043e-07

param="ESTIMATE"
d=read.table(paste0(holder,"_",param,"_summary"))
wilcox.test(d$cluster2,d$cluster4)
W = 1757, p-value = 0.03847
wilcox.test(d$cluster2,d$cluster1)
W = 5583, p-value = 2.851e-06
wilcox.test(d$cluster2,d$cluster3)
W = 1647, p-value = 0.0001015

param="IMPRES"
d=read.table(paste0(holder,"_",param,"_summary"))
wilcox.test(d$cluster2,d$cluster4)
W = 1722, p-value = 0.05412
wilcox.test(d$cluster2,d$cluster1)
W = 4707, p-value = 0.01871
wilcox.test(d$cluster2,d$cluster3)
W = 1446.5, p-value = 0.01422

param="IFNg"
d=read.table(paste0(holder,"_",param,"_summary"))
wilcox.test(d$cluster2,d$cluster4)
W = 1877, p-value = 0.004816
wilcox.test(d$cluster2,d$cluster1)
W = 5692, p-value = 6.354e-07
wilcox.test(d$cluster2,d$cluster3)
W = 1633, p-value = 0.0001551

param="Bindea"
d=read.table(paste0(holder,"_",param,"_summary"))
wilcox.test(d$cluster2,d$cluster4)
W = 1860, p-value = 0.006668
wilcox.test(d$cluster2,d$cluster1)
W = 5637, p-value = 1.37e-06
wilcox.test(d$cluster2,d$cluster3)
W = 1608, p-value = 0.0003227

param="age"
d=read.table(paste0(holder,"_",param,"_summary"))
wilcox.test(d$cluster2,d$cluster4)
W = 1860, p-value = 0.006668
wilcox.test(d$cluster2,d$cluster1)
W = 5637, p-value = 1.37e-06
wilcox.test(d$cluster2,d$cluster3)
W = 1608, p-value = 0.0003227


#perform correlation plot between #SNVs and immune scores
param="SNVs"
d=read.table(paste0(holder,"_",param,"_summary"))
param="ESTIMATE"
e=read.table(paste0(holder,"_",param,"_summary"))

cor.test(d$cluster1,e$cluster1)
-0.1220063, p-value = 0.1305
cor.test(d$cluster2,e$cluster2)
0.1527921, p-value = 0.2895
cor.test(d$cluster3,e$cluster3)
0.2163164, p-value = 0.1535
cor.test(d$cluster4,e$cluster4)
0.2060523, p-value = 0.1241


#nothing interesting....

########################
