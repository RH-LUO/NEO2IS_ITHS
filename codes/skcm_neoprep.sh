#!/bin/bash
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -p cn
indir=$HOME/lrh/SKCM/SNV/mutect2_results/filter.vcfs
if [ ! -e $HOME/lrh/SKCM/Neoantigen/Input.snp_vcfs/ ]; then 
mkdir $HOME/lrh/SKCM/Neoantigen  
mkdir $HOME/lrh/SKCM/Neoantigen/Input.indel_vcfs
mkdir $HOME/lrh/SKCM/Neoantigen/Input.snp_vcfs

cp $indir/*.indel.filtered.vcf $HOME/lrh/SKCM/Neoantigen/Input.indel_vcfs
cp $indir/*.snp.filtered.vcf $HOME/lrh/SKCM/Neoantigen/Input.snp_vcfs
#Rename
cd $HOME/lrh/SKCM/Neoantigen/Input.snp_vcfs
ls *.snp.filtered.vcf|while read id;do file=$(basename $id);sample=${file%%.*}; mv ${id} ${sample}.vcf;done                     

cd $HOME/lrh/SKCM/Neoantigen/Input.indel_vcfs
ls *.indel.filtered.vcf|while read id;do file=$(basename $id);sample=${file%%.*}; mv ${id} ${sample}.vcf; done                     
fi

source $HOME/miniconda3/bin/activate
conda activate python27 
outdir=$HOME/lrh/SKCM/Neoantigen
mkdir $outdir/SKCM.indel_results/
mkdir $outdir/SKCM.snp_results/
python $HOME/bin/NeoPredPipe/NeoPredPipe.py --preponly -I $HOME/lrh/SKCM/Neoantigen/Input.snp_vcfs -H $HOME/lrh/SKCM/HLA_I.txt -o $outdir/SKCM.snp_results/ -n SKCM -c 1 2 -E 8 9 10
python $HOME/bin/NeoPredPipe/NeoPredPipe.py --preponly -I $HOME/lrh/SKCM/Neoantigen/Input.indel_vcfs -H $HOME/lrh/SKCM/HLA_I.txt -o $outdir/SKCM.indel_results/ -n SKCM -c 1 2 -E 8 9 10
