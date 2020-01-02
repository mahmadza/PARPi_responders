
########################################
#compare parameters of the clusters
########################################

cd /people/mzabidi/tumor_project/data/ERpos_analysis


param=(age IMPRES IFNg Bindea ESTIMATE HRDscore SNVs indels dels ins)
colnum=(4 6 7 8 9 10 11 12 13 14)

x=$(mktemp)
for paramnum in ${!param[*]}; do
  for i in {1..4}; do
    cat summary_ward.D.txt | \
      awk -vy=${colnum[$paramnum]} -vp=cluster${i}.txt \
        'BEGIN{
          while((getline<p)>0) t[$1]=1
        }
        (NR>1 && $1 in t){print $1,$y}' > ${x}_cluster${i}_${param[$paramnum]}
  echo ${param[$paramnum]} $(awk '{print $2}' ${x}_cluster${i}_${param[$paramnum]} | grep -v NA | datamash median 1)

    done
done

#compile
for paramnum in ${!param[*]}; do
  maxrow=$(wc -l ${x}_cluster*_${param[$paramnum]} | sort -k1,1nr | awk 'NR==2{print $1}')
  for i in {1..4}; do
    cat ${x}_cluster${i}_${param[$paramnum]} | \
      awk -vp=$maxrow '{print $2} END{for(i=NR;i<p;i++) print "NA"}' > ${x}_cluster${i}_${param[$paramnum]}_OK
  done
  paste ${x}_cluster{1,2,3,4}_${param[$paramnum]}_OK | \
    awk -vOFS="\t" 'BEGIN{print "\tcluster1\tcluster2\tcluster3\tcluster4"}
      {print "dummy"++n,$0}' > ${x}_${param[$paramnum]}_summary
done


#calculate ANOVA here
for paramnum in ${!param[*]}; do
  pval=$(for i in {1..4}; do
  awk -vp=$i -vOFS="\t" '{print $2,p}' ${x}_cluster${i}_${param[$paramnum]}
  done | awk -vp=${param[$paramnum]} -vOFS="\t" 'BEGIN{print "\tage\tgroup"}{print NR,$0}' | \
  /people/mzabidi/tumor_project/utils/ANOVA_151019.R -i -)

  echo $paramnum ${param[$paramnum]} $pval
done
0 age 0.08612099
1 IMPRES 0.7130833
2 IFNg 0.2721241
3 Bindea 0.1611523
4 ESTIMATE 0.06838629
5 HRDscore 0.0003895215
6 SNVs 0.1519701
7 indels 0.0006812105
8 dels 2.816156e-05
9 ins 0.4948577


for paramnum in ${!param[*]}; do
  pval=$(for i in {1..4}; do
  awk -vp=$i -vOFS="\t" '{print $2,p}' ${x}_cluster${i}_${param[$paramnum]}
  done | awk -vp=${param[$paramnum]} -vOFS="\t" 'BEGIN{print "\tage\tgroup"}{print NR,$0}' | \
  /people/mzabidi/tumor_project/utils/ANOVA_151019.R -i -)

  echo ${param[$paramnum]} $pval
done
