#!/bin/bash
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -p cn
####CCF for SKCM###
source $HOME/miniconda3/bin/activate
data_path="$HOME/lrh/SKCM"
indir=$data_path/SNV/mutect2_results/filter.vcfs
# outdir=/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/CCF
outdir=/public/home/lorihan/lrh/Neoantigen_TCGA/SKCM/CCF.new
mkdir $outdir
# outfolder=$HOME/lrh/SKCM/ABSOLUTE
outfolder=$HOME/lrh/SKCM/ABSOLUTE.xena

ls $indir/vep.vcfs/*_vep.maf|while read id
do
file=$(basename $id )
sample=${file%%_*}
tumor=`echo ${sample:0:15}`
echo $id $sample $tumor
cat $id | sed '1,2d' | awk -F '\t' '{print $5"\t"$6"\t"$7"\t"$1":"$16"\t"$37"\t"$46"\t"$41"\t"$42"\t"2"\t"0"\t"2"\t"}' >$outdir/${tumor}.tmp.tsv
done

# data_path="$HOME/lrh/Neoantigen_TCGA/SKCM"
data_path="$HOME/lrh/SKCM"

cat $data_path/TCGA_samples_maf_cn.txt|while read sample
do
echo $sample
# cat $outfolder/${sample}.seg | sed '1d' | awk 'BEGIN{OFS="\t"}{print $0"\t"((2^$6)*2)}'| awk 'BEGIN{OFS="\t"}{if ($7!=0)print $0}' | cut -f 2-7  >$outdir/${sample}.bed
cat $outfolder/${sample}.seg | sed '1d' | awk 'BEGIN{OFS="\t"}{print $0"\t"((2^$5)*2)}'| awk 'BEGIN{OFS="\t"}{if ($6!=0)print $0}' | cut -f 2-6  >$outdir/${sample}.bed
bedtools window -a $outdir/${sample}.tmp.tsv -b $outdir/${sample}.bed | cut -f 4-10,16 | awk 'BEGIN{OFS="\t";print "mutation_id\tHGVSp_Short\tall_effects\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn"}{print $0}' >$outdir/${sample}.tsv
done
rm $outdir/*.tmp.tsv





# 
# ####CCF for TCGA-LUNG###
# data_path="$HOME/lrh/Neoantigen_NSCLC"
# outdir=$data_path/TCGA
# #cp $HOME/lrh/LUAD/SNV/mutect2_results/filter.vcfs/vep.vcfs/*_vep.maf $outdir
# #cp $HOME/lrh/LUSC/SNV/mutect2_results/filter.vcfs/vep.vcfs/*_vep.maf $outdir
# 
# ls $outdir/*_vep.maf|while read id
# do
# file=$(basename $id )
# sample=${file%%_*}
# echo $id $sample
# cat $id | sed '1,2d' | awk -F '\t' '{print $5"\t"$6"\t"$7"\t"$1":"$16"\t"$37"\t"$46"\t"$41"\t"$42"\t"2"\t"0"\t"2"\t"}' >$outdir/${sample}.tmp.tsv
# done
# 
# data_path="$HOME/lrh/Neoantigen_NSCLC"
# outdir=$data_path/TCGA
# 
# cat $data_path/TCGA_samples_maf_cn.txt|while read sample
# do
# echo $sample
# cat $outdir/${sample}.seg | sed '1d' | awk 'BEGIN{OFS="\t"}{print $0"\t"((2^$6)*2)}'| awk 'BEGIN{OFS="\t"}{if ($7!=0)print $0}' | cut -f 2-7  >$outdir/${sample}.bed
# bedtools window -a $outdir/${sample}.tmp.tsv -b $outdir/${sample}.bed | cut -f 4-10,17 | awk 'BEGIN{OFS="\t";print "mutation_id\tHGVSp_Short\tall_effects\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn"}{print $0}' >$outdir/${sample}.tsv
# done
# rm *.tmp.tsv