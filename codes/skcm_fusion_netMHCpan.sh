#!/bin/sh
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -p cn
#source $HOME/miniconda3/bin/activate
infolder1=$HOME/lrh/SKCM/Neoantigen/Input.fusion_tsvs
outfolder=$HOME/lrh/SKCM/Neoantigen/SKCM.fusion_results
ls $infolder1/*.fasta|while read id
do
file=$(basename $id )
name=${file%%.*}
echo $name
hla_line=`cat $HOME/lrh/SKCM/Neoantigen/SKCM.fusion_results/hla/${name}.hla`
if [ ! -e $outfolder/${name}.fusion.sb ]; then 
echo $hla_line
echo -e "\n*** SKCM_fusion_netMHCpan started at $(date +'%T %F') ***\n"
netMHCpan -BA -s -a $hla_line $infolder1/${name}.fasta>$outfolder/${name}.fusion.out
less -S $outfolder/${name}.fusion.out|grep "Pos"|head -1 >$outfolder/header
less -S $outfolder/${name}.fusion.out|grep "<= SB"|awk '{gsub("<= ","<=");print $0}'>$outfolder/${name}.fusion.pep
cat $outfolder/header $outfolder/${name}.fusion.pep >$outfolder/${name}.fusion.sb
fi
done
