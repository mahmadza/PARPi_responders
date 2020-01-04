#######################################################################################
#investigate the differential genes expressed in the different clusters
#use pamr R package to explore whether gene expression can separate these clusters
#######################################################################################

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







########################
