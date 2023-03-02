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
mkdir $outdir1/sample.id
mkdir $outdir/vep.vcfs
GATK_bundle="$HOME/lrh/hg38"
ref="$GATK_bundle/Homo_sapiens_assembly38.fasta"

cat $data_path/TCGA_SKCM_VCF_MuTect2.txt|awk 'NR==1{next}{print $0}'|while read id
do
folder=`echo $id | awk '{print $3 }'`
file=`echo $id | awk '{print $4}'`
tumor=`echo $id | awk '{print $1}'`
normal=`echo $id | awk '{print $5}'`

bcftools query -l $outdir1/mutect2_results/${tumor}.vcf>$outdir1/sample.id/${tumor}_id
cat $outdir1/sample.id/${tumor}_id|sed -n '1p'>$outdir1/sample.id/${normal}_sample_name
cat $outdir1/sample.id/${tumor}_id|sed -n '2p'>$outdir1/sample.id/${tumor}_sample_name

tumor_SM=`cat $outdir1/sample.id/${tumor}_sample_name`
normal_SM=`cat $outdir1/sample.id/${normal}_sample_name`

echo $name $tumor_SM $normal_SM    
echo "start vep_annotation for ${tumor} " `date`

if [ ! -e $outdir/vep.vcfs/${tumor}_vep.vcf ]; then   
echo $name $sample $tumor
echo "start vep_annotation for ${id} " `date`
vep --cache --offline --format vcf --vcf --force_overwrite --everything  --fork 4 --assembly GRCh38 \
 --dir_cache $HOME/lrh/vep_data/ \
 --dir_plugins $HOME/lrh/vep_data/Plugins/ \
 --input_file $outdir/${tumor}.combine.filtered.vcf \
 --output_file $outdir/vep.vcfs/${tumor}_vep.vcf \
 --fasta $ref \
 --plugin Downstream \
 --plugin Wildtype \
 --symbol \
 --terms SO \
 --tsl 
echo "end vep_annotation for ${id} " `date`

echo "start vep2maf for ${id} " `date`
echo $name $tumor_SM $normal_SM     
         perl $HOME/lrh/vep_data/vcf2maf.pl \
        --input-vcf $outdir/vep.vcfs/${tumor}_vep.vcf \
        --output-maf $outdir/vep.vcfs/${tumor}_vep.maf \
        --ref-fasta $ref \
        --tumor-id  ${tumor_SM} \
        --normal-id ${normal_SM} \
        --ncbi-build GRCh38 \
        --tmp-dir //tmp \
        --inhibit-vep
echo "end vep2maf for ${id} " `date`
fi
done

outdir=$data_path/SNV/mutect2_results/filter.vcfs
if [ `ls $outdir/vep.vcfs/*vep.maf|wc -l`==`ls $outdir1/mutect2_results/*.vcf|wc -l` ]; then
cat $outdir/vep.vcfs/*maf | grep -v '^#'| grep -v '^Hugo_Symbol' >$outdir/vep.vcfs/tmp 
id=`cat $data_path/TCGA_SKCM_VCF_MuTect2.txt|awk 'NR==1{next}{print $0}'|head -1`
tumor=`echo $id | awk '{print $1}'`
grep 'Hugo_Symbol' $outdir/vep.vcfs/${tumor}_vep.maf >$outdir/vep.vcfs/header
cat $outdir/vep.vcfs/header $outdir/vep.vcfs/tmp >$outdir/vep.vcfs/vep_merge.maf
cd $HOME/lrh/SKCM/SNV/mutect2_results/filter.vcfs/vep.vcfs
cut -f 5-7,12,13,1,16 vep_merge.maf |cut -f 2-7  > 1
cut -f 5-7,12,13,1,16 vep_merge.maf |cut -f 1 > 2
paste 1 2 > for_annovar.input 

tool_path="$HOME/lrh/bin"
data_path="$HOME/lrh/SKCM"
bin=$HOME/lrh/bin/annovar
db=$HOME/lrh/bin/annovar/humandb/
indir=$HOME/lrh/SKCM/SNV/mutect2_results/filter.vcfs/vep.vcfs
echo "start annoFilter for SKCM " `date`
perl $bin/table_annovar.pl $indir/for_annovar.input $db -buildver hg38 -out $indir/tmp -protocol genomicSuperDups,rmsk,1000g2015aug_all,exac03,gnomad_genome -operation r,r,f,f,f -nastring NA
echo "annoFilter end" `date`
fi

