##############################################################
#      cluster ER+ samples using mutational signatures
##############################################################

cd ~/work/tumor_project/data/redo_samples/9_Strelka2_indels

#grab ER+ samples
x=$(mktemp)
cat ../../../annotations/Study_metadata_8.5.19.csv | \
  awk -vFS="," '{print $1,$8}' | \
    awk '$2=="P"{print $1}' > ${x}_ERpos
#352 samples

#grab signatures, use it for clustering
( l=SD0012
  cat deconstructSigs_signatures_breast_only/${l}.txt | \
  awk 'NR==1{line=$0}
      NR==3 || NR==5{
        for(i=1;i<=NF;i++)
          line=line"\t"$i
      }
      END{print line}'
  cat ${x}_ERpos | while read l; do
    if [ -s deconstructSigs_signatures_breast_only/${l}.txt ]; then
      cat deconstructSigs_signatures_breast_only/${l}.txt | \
        awk 'NR==2{line=$0}
            NR==4 || NR==6{
              for(i=2;i<=NF;i++)
                line=line"\t"$i
            }
            END{print line}'
    fi
  done ) > ${x}_heatmap


#enter R here, or copy and paste these codes to another window with R running
R

library(gplots)
paste0 <- function( ..., sep="" ) paste( ..., sep = sep )
holder="/tmp/tmp.gogJCAUBoA"

#n breaks
#n-1 colors
my.breaks=seq(0,0.8,length.out=101)       #101 breaks
my.colors<-colorpanel(100,rgb(49,105,154,25,maxColorValue=255),
rgb(255,255,255,25,maxColorValue=255),rgb(164,25,32,25,maxColorValue=255))       #100 colors

e<-read.table(paste0(holder,"_heatmap"))
my.cex=0.5
for(curr_method in c("ward.D","ward.D2"))
{
  pdf(paste0("/home/mzabidi/Erpos_samples_",curr_method,".pdf"))
  hclust2 <- function(x, method=curr_method, ...)
    hclust(x, method=method, ...)
  title=paste0("ERpos sample signatures\nmethod = ",curr_method)
  v=heatmap.2(as.matrix(t(e)),main=title,dendrogram='both',trace="none",
  cexCol=0.08,cexRow=0.8,
    density.info="none",key=TRUE,symkey=TRUE,
    breaks=my.breaks,scale="none",col=my.colors,
    hclustfun=hclust2)

  dev.off()

  #print out the sample order
  fileConn<-file(paste0(holder,"_order_",curr_method))
  writeLines(c(row.names(v$carpet)), fileConn)
  close(fileConn)

}


#


###################################
#    plot HRD score heatmap
###################################

method=ward.D
( echo -en "\tHRDscore\tHRDscore_dummy\n"
  cat ${x}_order_${method} | while read l; do

    HRD=$(cat ../../../annotations/HRDscores.txt | awk -vp=$l 'BEGIN{c="-2"}($1==p){c=$5} END{print c}' )

    echo -en $l"\t"$HRD"\t-2\n"

  done ) > ${x}_HRD_${method}


#enter R here, or copy and paste these codes to another window with R running
R
my.cex=0.5
curr_method="ward.D"
e<-read.table(paste0(holder,"_HRD_",curr_method))
my.breaks=c(min(e),seq(0,100,length.out=10))
my.colors<-colorpanel(10,"lavenderblush","darkred")

pdf(paste0("/home/mzabidi/R_output/HRDscores_",curr_method,".pdf"))
heatmap.2(as.matrix(t(e)),main=title,trace="none",
cexCol=0.3,cexRow=0.8,
  density.info="none",key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors,
            dendrogram='none',Rowv=FALSE,Colv=FALSE)
dev.off()


###################################
#    plot IntClust as heatmap
###################################
( echo -en "\tIntClust\tIntClust_dummy\n"
  cat ${x}_order_${method} | while read l; do

    IntClust=$(cat ../../../annotations/Study_metadata_8.5.19.csv | \
      awk -vFS="," -vp=$l 'BEGIN{c="-2"}($1==p){c=$18} END{print c}' | \
        awk '$1=="4-"{$1="4.3"}$1=="4+"{$1="4.8"}$1=="NA"{$1="-2"}{print}')
    #IntClust4- -> 4.3
    #IntClust4+ -> 4.8
    echo -en $l"\t"$IntClust"\t-2\n"

  done
  #also add extra columns as calibrator
  for i in -2 {1..3} 4.3 4.8 {5..10}; do
    echo -en "IntClust"$i"\t"$i"\t-2\n"
  done
  ) > ${x}_IntClust_${method}


#enter R here, or copy and paste these codes to another window with R running
R
curr_method="ward.D"
e<-read.table(paste0(holder,"_IntClust_",curr_method))
library(RColorBrewer)
my.breaks=c(-2,0,1,2,3,4.3,4.8,5,6,7,8,9,10)       #12 breaks
my.colors=c("grey",brewer.pal(11,"Paired"))       #"

pdf(paste0("/home/mzabidi/R_output/IntClust_",curr_method,".pdf"))

heatmap.2(as.matrix(t(e)),main=title,trace="none",
cexCol=0.3,cexRow=0.8,
  density.info="none",key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors,
            dendrogram='none',Rowv=FALSE,Colv=FALSE)
legend("bottomleft",c("NA",paste0("IntClust",c(1,2,3,"4-","4+",5,6,7,8,9,10))),col=my.colors,
      bty="n",pch=15)
#"
dev.off()


###################################
#    plot PAM50 as heatmap
###################################
( echo -en "\PAM50\tPAM50_dummy\n"
  cat ${x}_order_${method} | while read l; do

    pam50=$(cat ../../../annotations/Study_metadata_8.5.19.csv | \
      awk -vFS="," -vp=$l 'BEGIN{c="-2"}($1==p){c=$17} END{print c}' | \
        awk '$1=="Normal"{c=1}
              $1=="LumA"{c=2}
              $1=="LumB"{c=3}
              $1=="Her2"{c=4}
              $1=="Basal"{c=5}
              $1=="NA"{c="-2"}
              END{print c}')
    echo -en $l"\t"$pam50"\t-2\n"
  done

  #also add extra columns as calibrator
  echo -en "Normal\t1\t-2\n"
  echo -en "LumA\t2\t-2\n"
  echo -en "LumB\t3\t-2\n"
  echo -en "Her2\t4\t-2\n"
  echo -en "Basal\t5\t-2\n"
  echo -en "Class_NA\t-2\t-2\n"
  ) > ${x}_PAM50_${method}

#enter R here, or copy and paste these codes to another window with R running
R
curr_method="ward.D"
e<-read.table(paste0(holder,"_PAM50_",curr_method))
library(RColorBrewer)
my.breaks=c(-2,0,1,2,3,4,5)       #6 breaks
my.colors=c("grey",brewer.pal(5,"Set2"))       #"5 breaks

pdf(paste0("/home/mzabidi/R_output/PAM50_",curr_method,".pdf"))
heatmap.2(as.matrix(t(e)),main=paste0("PAM50 clustering ",curr_method),trace="none",
cexCol=0.3,cexRow=0.8,
  density.info="none",key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors,
            dendrogram='none',Rowv=FALSE,Colv=FALSE)
#"
legend("bottomleft",c("NA","Normal","Luminal A","Luminal B","Her2-enriched","Basal"),col=my.colors,
      bty="n",pch=15)
dev.off()


#"
################################################
#    plot BRCA1/2 germline status as heatmap
################################################
( echo -en "\germline_status\tgermline_status_dummy\n"
  cat ${x}_order_${method} | while read l; do
    carrier=$(cat ../../../annotations/All_HR_carriers.txt | \
      awk -vp=$l 'BEGIN{c="NA"}$1==p{c=$2} END{print c}' | \
        awk 'BEGIN{d="-2"}
        $1~/^BRCA1_del/{d=1}
            $1~/^BRCA2_del/{d=2}
            $1~/^PALB2_del/{d=3}
            $1~/^ATMdel/{d=4}
            $1~/^CHEK2_del/{d=5}
            END{print d}')
    echo -en $l"\t"$carrier"\t-2\n"
  done

  #also add extra columns as calibrator
  echo -en "BRCA1_del\t1\t-2\n"
  echo -en "BRCA2_del\t2\t-2\n"
  echo -en "PALB2_del\t3\t-2\n"
  echo -en "ATM_del\t4\t-2\n"
  echo -en "CHEK2_del\t5\t-2\n"
  echo -en "Class_NA\t-2\t-2\n"

  ) > ${x}_germssss_${method}



curr_method="ward.D"
e<-read.table(paste0(holder,"_germssss_",curr_method))
library(RColorBrewer)
my.breaks=c(-2,0,1,2,3,4,5)       #6 breaks
my.colors=c("grey",brewer.pal(5,"Set2"))       #"5 breaks

pdf(paste0("/home/mzabidi/R_output/carrierstatus_",curr_method,".pdf"))
heatmap.2(as.matrix(t(e)),main=paste0("germline status clustering ",curr_method),trace="none",
cexCol=0.3,cexRow=0.8,
  density.info="none",key=TRUE,symkey=TRUE,
  breaks=my.breaks,scale="none",col=my.colors,
            dendrogram='none',Rowv=FALSE,Colv=FALSE)
legend("bottomleft",c("NA","BRCA1","BRCA2","PALB2","ATM","CHEK2"),col=my.colors,
      bty="n",pch=15)

dev.off()



#"
#do this in Shell...
( echo -en "\tPAM50\tcarrier\tIntClust\tIMPRES\tHRDscore\tSNVs\tindels\tdels\tins\n"
  cat ${x}_order_${method} | while read l; do

    carrier=$(cat ../../../annotations/All_HR_carriers.txt | \
      awk -vp=$l 'BEGIN{c="NA"}$1==p{c=$2} END{print c}' | \
        awk 'BEGIN{d="NA"}
        $1~/^BRCA1_del/{d=5}
            $1~/^BRCA2_del/{d=4}
            $1~/^PALB2_del/{d=3}
            $1~/^ATMdel/{d=2}
            $1~/^CHEK2_del/{d=1}
            END{print d}')
    PAM50=$(cat ../../../annotations/immune_score/SJMC_immunescores.txt | \
      awk -vp=$l 'BEGIN{c="NA"}($1==p){c=$3} END{print c}' )
    IntClust=$(cat ../../../annotations/immune_score/SJMC_immunescores.txt | \
      awk -vp=$l 'BEGIN{c="NA"}($1==p){c=$4} END{print c}' )
    IMPRES=$(cat ../../../annotations/immune_score/SJMC_immunescores.txt | \
      awk -vp=$l 'BEGIN{c="NA"}($1==p){c=$5} END{print c}' )
    HRD=$(cat ../../../annotations/HRDscores.txt | awk -vp=$l 'BEGIN{c="NA"}($1==p){c=$5} END{print c}' )

    if [ -s ${l}.filtered.vcf.gz ]; then
      SNV=$(gunzip -c ${l}.filtered.vcf.gz | awk '$0!~/^#/' | awk 'length($4)==length($5){c++} END{print c+0}')
      indels=$(gunzip -c ${l}.filtered.vcf.gz | awk '$0!~/^#/' | awk 'length($4)!=length($5){c++} END{print c+0}')
      dels=$(gunzip -c ${l}.filtered.vcf.gz | awk '$0!~/^#/' | awk 'length($4)>length($5){c++} END{print c+0}')
      ins=$(gunzip -c ${l}.filtered.vcf.gz | awk '$0!~/^#/' | awk 'length($4)<length($5){c++} END{print c+0}')
    else
      SNV="NA"
      indels="NA"
      dels="NA"
      ins="NA"
    fi

  echo -en $l"\t"$PAM50"\t"$carrier"\t"$IntClust"\t"$IMPRES"\t"$HRD"\t"$SNV"\t"$indels"\t"$dels"\t"$ins"\n"

done ) > ${x}_summary_${method}



#plot as barplots
d<-read.table(paste0(holder,"_summary_",curr_method))
pdf(paste0("/home/mzabidi/R_output/ERpos_samples_",curr_method,"_panels.pdf"))

par(mfrow=c(4,2))
for(i in seq(2,9))
{
  barplot(d[,i],names.arg=row.names(d),main=colnames(d)[i],
  col=rainbow(9)[i-1],border=rainbow(9)[i-1],cex.names=0.8,las=2)
}

dev.off()



########################
