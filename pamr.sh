




#investigate the differential genes expressed in the different clusters
#use pamr R package to explore whether gene expression can separate these clusters

cd ~/work/tumor_project/data/ERpos_analysis

#from the heatmap, manually determine the cluster boundaries
#and appropriately separate the samples into clusters
#Cluster1: high signature 1: SD1442 till SD0029
#Cluster2: high signatures 2&13: SD0899 till SD0164
#cluster3: mixed: SD1009 till SD0171
#cluster4: high signatures 3: SD0734 till SD1321

cat /tmp/tmp.gogJCAUBoA_order_ward.D | \
  awk '(!n){print}($1=="SD0029"){n=1}' > cluster1.txt
cat /tmp/tmp.gogJCAUBoA_order_ward.D | \
  awk 'BEGIN{while((getline<"cluster1.txt")>0) k[$1]=1}
    !($1 in k)' | \
    awk '(!n){print}($1=="SD0164"){n=1}' > cluster2.txt
cat /tmp/tmp.gogJCAUBoA_order_ward.D | \
  awk 'BEGIN{while((getline<"cluster1.txt")>0) k[$1]=1}
    !($1 in k)' | \
    awk 'BEGIN{while((getline<"cluster2.txt")>0) k[$1]=1}
      !($1 in k)' | \
          awk '(!n){print}($1=="SD0171"){n=1}' > cluster3.txt
cat /tmp/tmp.gogJCAUBoA_order_ward.D | \
  awk 'BEGIN{while((getline<"cluster1.txt")>0) k[$1]=1}
    !($1 in k)' | \
    awk 'BEGIN{while((getline<"cluster2.txt")>0) k[$1]=1}
      !($1 in k)' | \
      awk 'BEGIN{while((getline<"cluster3.txt")>0) k[$1]=1}
        !($1 in k)' > cluster4.txt

159 cluster1.txt
 55 cluster2.txt
 49 cluster3.txt
 62 cluster4.txt


#grab the gene expression from each sample
rna=$(mktemp)
for num in {1..4}; do
  cat cluster${num}.txt | while read l; do
    echo $l
    #grab column number of the sample, if exist
    #if doesn't exist, skip
    colnum=$(zcat ../../annotations/SD_RSEM_Output_TPM_3.txt.gz | \
      awk -vp=$l 'NR==1{
        for(i=4;i<=NF;i++)
          if($i==p)
          {
            print i
            exit
          }
          print "NA"      #if RNA-seq not found
      }')

      #grab gene expression if exists
      if [ "$colnum" != "NA" ]; then
        zcat ../../annotations/SD_RSEM_Output_TPM_3.txt.gz | \
          cut -f $colnum > ${rna}_cluster${num}_${l}
      fi
  done &
done

#combine samples from the same clusters
#and create headers (i.e cluster numbers)
for num in {1..4}; do
  paste ${rna}_cluster${num}_SD* > ${rna}_cluster${num}
  paste ${rna}_cluster${num}_* | \
    awk -vp=$num 'NR==1{
      line=p
      for(i=2;i<=NF;i++)
        line=line"\t"p
      print line
    }' > ${rna}_cluster${num}_header
done

#create gene expression matrix
paste <(zcat ../../annotations/SD_RSEM_Output_TPM_3.txt.gz | \
  awk 'NR==1{print ""}NR>1{print $2}') ${rna}_cluster{1..4} > ${rna}_cluster_complete
#57906 genes

#normalize by z-score
#first, find mean and stdev
awk 'NR>1' ${rna}_cluster_complete | \
  awk -vOFS="\t" 'NR>1{
    sum=0
    sumsq=0
    count=0
    zero_count=0
    max=$2
    min=$2
    for(i=2;i<=NF;i++)
      if($i!="NA")
      {
        sum+=$i
        sumsq+=$i*$i
        count++
        if($i>max) max=$i
        if($i<min) min=$i
        if($i=="0") zero_count++
      }
    mean=sum/count
    stdev=sqrt(sumsq/count-(sum/count)^2)
    if(mean)
      print $1,sum,mean,stdev,stdev/mean*100,min,max,zero_count
    else
      print $1,sum,mean,stdev,"NA",min,max,zero_count
  }' > ${rna}_raw_summary
#gene_name,sum,mean,stdev,%stdev,min,max,#zeros

#compute z-transformed gene expression
#z=(x-mean)/stdev
#exclude those with stdev=0 (actually none is excluded by this)
#and mean tpkm<=2
#also print gene names
awk 'NR>1' ${rna}_cluster_complete | \
  awk -vp=${rna}_raw_summary -vq=${rna}_cluster_genenames \
    'BEGIN{
      while((getline<p)>0)
      {
        mean[$1]=$3
        stdev[$1]=$4
        zero_count[$1]=$8
      }
    }
    (stdev[$1] && mean[$1]>=2){
      line=($2-mean[$1])/stdev[$1]
      for(i=3;i<=NF;i++)
        line=line"\t"($i-mean[$1])/stdev[$1]
      print line
      print $1 > q
    }' > ${rna}_cluster_complete_zcore_trans
#14639 genes

paste ${rna}_cluster{1..4}_header > ${rna}_cluster_header_summary

echo $rna
/tmp/tmp.laukMitnA2




#run pamr package to find gene signatures separating the 4 clusters
R

#load gplots for drawing heatmap
#load pamr for analysis
library(gplots)
library(pamr)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )
holder="/tmp/tmp.b4keRzvD7u"

#load z-transform gene expression into matrix
exp=as.matrix(read.delim(paste0(holder,"_cluster_complete_zcore_trans"),header=F))
#load headers (i.e cluster numbers)
header=as.integer(read.table(paste0(holder,"_cluster_header_summary")))
#load gene names
gN=as.character(read.table(paste0(holder,"_cluster_genenames"))$V1)
#""
#input data ready for classification
#also input gene names and ID
mydata=list(x=exp,y=factor(header),genenames=gN,geneid=gN)

#""
#train classifier
mytrain=pamr.train(mydata)
#cross validate the classifier at 10-fold
myresults=pamr.cv(mytrain,mydata)
#plot the cross-validated error curves
pdf(paste0("/home/mzabidi/R_output/missclassification_errors.pdf"))
pamr.plotcv(myresults)
dev.off()
#""
#here, it shows that
#1)there'll be little variation in classification
#generally variation in classification will be lower if higher threshold
#2)for the different clusters,
#the misclassification error will increase as increase the threshold values
#EXCEPT for Cluster 1
#cluster 4 has high errors no matter what

#compute the confusion matrix for a particular model
thold=1.15   #gives the best cross-validated class probabilities plot
#plot the cross-validated class probabilities by class
pdf(paste0("/home/mzabidi/R_output/crossvalidated_prob_",thold,".pdf"))
pamr.plotcvprob(myresults,mydata,threshold=thold)
dev.off()
#it seems that we cant get cluster 4



##""
#plot the class centroids
pdf(paste0("/home/mzabidi/R_output/classcentroids_",thold,".pdf"))
pamr.plotcen(mytrain,mydata,threshold=thold)
dev.off()
#just means how discriminatory each gene is
#cluster4 in general is not looking really good


#plot of the most significant genes
#change the threshold since otherwise will plot too many genes
#save into pdfs
#increase thold so that can at least plot some genes
thold_b=2.75
pdf(paste0("/home/mzabidi/R_output/pamr_genes_",thold_b,".pdf"))
pamr.geneplot(mytrain,mydata,threshold=thold_b)
dev.off()


#list significant genes
genes=pamr.listgenes(mytrain,mydata,threshold=thold,genenames=T)
head(genes,20)
id       name     1-score   2-score  3-score 4-score
[1,] "CTSW"   "CTSW"   "-0.0413" "0.2404" "0"     "0"
[2,] "ASRGL1" "ASRGL1" "-0.0508" "0.2366" "0"     "0"
[3,] "IL18BP" "IL18BP" "-0.0544" "0.2312" "0"     "0"
[4,] "GBP1"   "GBP1"   "-0.0454" "0.2277" "0"     "0"
[5,] "RAB20"  "RAB20"  "-0.0517" "0.2245" "0"     "0"
[6,] "CLCA2"  "CLCA2"  "-0.021"  "0.2234" "0"     "0"
[7,] "GBP2"   "GBP2"   "-0.0329" "0.2202" "0"     "0"
[8,] "WARS"   "WARS"   "-0.0496" "0.2066" "0"     "0"
[9,] "DBI"    "DBI"    "-0.0013" "0.2035" "0"     "0"
[10,] "ZNF296" "ZNF296" "-0.003"  "0.2013" "0"     "-0.0067"
[11,] "CXCL9"  "CXCL9"  "-0.0401" "0.198"  "0"     "0"
[12,] "GBP5"   "GBP5"   "-0.0508" "0.1969" "0"     "0"
[13,] "CASP4"  "CASP4"  "-0.0558" "0.1958" "0"     "0"
[14,] "CD2"    "CD2"    "-0.0414" "0.1933" "0"     "0"
[15,] "SPTLC2" "SPTLC2" "-0.0255" "0.1901" "0"     "0"
[16,] "IRF1"   "IRF1"   "-0.0425" "0.1895" "0"     "0"
[17,] "PTPN7"  "PTPN7"  "-0.0406" "0.1851" "0"     "0"
[18,] "CD8A"   "CD8A"   "-0.055"  "0.1851" "0"     "0"
[19,] "EAF2"   "EAF2"   "0"       "0.1844" "0"     "0"
[20,] "RAB27A" "RAB27A" "-0.0147" "0.184"  "0"     "0"








########################
