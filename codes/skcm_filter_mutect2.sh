#!/bin/bash
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -p cn
source $HOME/miniconda3/bin/activate
tool_path="$HOME/lrh/bin"
data_path="$HOME/lrh/SKCM"
indir=$data_path/skcm_mutect2
outdir=$data_path/SNV/mutect2_results/filter.vcfs
outdir1=$data_path/SNV

GATK_bundle="$HOME/lrh/hg38"
ref="$GATK_bundle/Homo_sapiens_assembly38.fasta"

cat $data_path/TCGA_SKCM_VCF_MuTect2.txt|awk 'NR==1{next}{print $0}'|while read id
do
folder=`echo $id | awk '{print $3 }'`
file=`echo $id | awk '{print $4}'`
tumor=`echo $id | awk '{print $1}'`

echo $folder $file $tumor
cp $indir/$folder/$file $outdir1/mutect2_results/${tumor}.vcf.gz
done

gunzip $outdir1/mutect2_results/*.gz
ls $outdir1/mutect2_results/*.gz|while read id
do
file=`echo $id | awk '{gsub(".gz","");print $0}'`
echo $id $file
mv -i $id $file
done

cat $data_path/TCGA_SKCM_VCF_MuTect2.txt|awk 'NR==1{next}{print $0}'|while read id
do
folder=`echo $id | awk '{print $3 }'`
file=`echo $id | awk '{print $4}'`
tumor=`echo $id | awk '{print $1}'`

if [ ! -e $outdir/${tumor}.snp.vcf ]; then   
echo "start vcf_filter for ${tumor} " `date`
##SNP Filter
gatk SelectVariants -select-type SNP -V $outdir1/mutect2_results/${tumor}.vcf -O $outdir/${tumor}.snp.vcf

gatk VariantFiltration -V $outdir/${tumor}.snp.vcf --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || QUAL < 30.0" --filter-name "Filter" -O  $outdir/${tumor}.snp.filter.vcf

##
#`--exclude-filtered true` filtering with pass
gatk SelectVariants --exclude-filtered true -V  $outdir/${tumor}.snp.filter.vcf -O  $outdir/${tumor}.snp.filtered.vcf

##INDEL Filter
gatk SelectVariants -select-type INDEL -V $outdir1/mutect2_results/${tumor}.vcf -O  $outdir/${tumor}.indel.vcf

gatk VariantFiltration -V  $outdir/${tumor}.indel.vcf --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || QUAL < 30.0 || ReadPosRankSum < -20.0" --filter-name "Filter" -O $outdir/${tumor}.indel.filter.vcf


#`--exclude-filtered true` filtering with pass
gatk SelectVariants --exclude-filtered true -V $outdir/${tumor}.indel.filter.vcf -O $outdir/${tumor}.indel.filtered.vcf
##
gatk MergeVcfs -I $outdir/${tumor}.snp.filtered.vcf -I $outdir/${tumor}.indel.filtered.vcf -O $outdir/${tumor}.combine.filtered.vcf

gatk MergeVcfs -I $outdir/${tumor}.snp.filter.vcf -I $outdir/${tumor}.indel.filter.vcf -O $outdir/${tumor}.combine.filter.vcf
fi
done
