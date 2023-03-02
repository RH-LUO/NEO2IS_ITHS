#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p cn
source $HOME/miniconda3/bin/activate
indir=$HOME/lrh/SKCM/Neoantigen/SKCM.snp_results
outdir=$HOME/lrh/SKCM/Neoantigen/SKCM.snp_results/fastaFiles

rm $indir/fasta.missed.file
cd $outdir
ls $outdir/*.reformat.fasta|while read id
do
file=$(basename $id )
tumor=${file%%.*}
echo $file $tumor

checkfilesize()
{

       time=$(date "+%Y-%m-%d %H:%M:%S")   #system time
       filename=${tumor}.reformat.fasta
       echo $id $filename
       if [ ! -f $filename ];then          #make a file
             # touch $outdir/$filename
              touch $indir/fasta.missed.file
             # echo "" >> $outdir/$filename
          #    echo "" >> ../fasta.log
          fi
          filezise=`ls -l $filename|awk '{print $5}'`
          maxsize=$((0))  
          echo $filesize $maxsize                      #max memory 10k
         if [ ! -s $filename ];then   
         echo "$tumor" >> $indir/fasta.missed.file #make a copy then delete
         rm $filename
           fi
}
checkfilesize
done

rm -r $outdir/test.fa
mkdir $outdir/test.fa

ls $indir/fastaFiles/*.reformat.fasta|while read id
do
file=$(basename $id )
name=${file%%.*}
echo $name
cat $indir/fastaFiles/${name}.reformat.fasta|awk NF|perl -p -i -e "s/\n//g" > $outdir/test.fa/${name}.test.fa
done

rm -r $indir/test/
mkdir $indir/test/

mkdir $indir/test/test.fasta
mkdir $indir/test/pos
mkdir $indir/test/netMHCpan

outdir=$HOME/lrh/SKCM/Neoantigen/SKCM.snp_results/test/test.fasta
cd $outdir

echo -e "\n*** SKCM_neopeptide_filter.R started at $(date +'%T %F') ***\n"
Rscript $HOME/R/rscripts/SKCM_neopeptide_filter.R

rm $indir/test/fasta.file
ls $outdir/*.test.fasta |while read id 
do
checkfilesize()
{
         file=$(basename $id )
         name=${file%%.*}
         filename=${name}.test.fasta
    echo $name
         touch $indir/test/fasta.file
          filezise=`ls -l $filename|awk '{print $5}'`
          maxsize=$((0))  
          echo $filesize $maxsize                      
          if [ ! -s $filename ];then   
         echo "$name" >> $indir/test/fasta.file #
         rm  $filename
           fi
}
checkfilesize
done

ls $indir/test/test.fasta/*.fasta|while read id
do
file=$(basename $id )
name=${file%%.*}
echo $name
cat $indir/test/test.fasta/${name}.test.fasta|grep -v "silent"|awk '{gsub(")",")\n"); print $0}'>$indir/test/test.fasta/${name}.snp.fa
  
cat  $indir/test/test.fasta/${name}.snp.fa|grep ">" >$indir/test/pos/${name}.snp.sites
cat $indir/test/pos/${name}.snp.sites| awk  'BEGIN { FS = "position" } ; { print $2 }'|awk  'BEGIN { FS = "changed" } ; { print $1 }'|awk 'BEGIN{OFS=FS="\t"}{pos=$1;pos_up=pos+10; pos_dw=pos-10;if(pos_dw<0) pos_dw=0;print pos, pos_dw, pos_up;}'  >$indir/test/pos/${name}.snp.range

cat  $indir/test/pos/${name}.snp.sites|awk -F ' ' '{print $1}' |awk '{gsub(">","");print $0}'>$indir/test/pos/${name}.gene
paste $indir/test/pos/${name}.gene  $indir/test/pos/${name}.snp.range > $indir/test/pos/${name}.snp.bed
done

echo -e "\n*** SKCM_neopeptide_shorten.R started at $(date +'%T %F') ***\n"
Rscript $HOME/R/rscripts/SKCM_neopeptide_shorten.R

mkdir $indir/test/hla
Rscript $HOME/R/rscripts/skcm_hla.R


ls $indir/test/test.fasta/*.snp.fa|while read id
do
file=$(basename $id )
name=${file%%.*}
echo $name
#hla_line=`cat $id`
hla_line=`cat $indir/test/hla/${name}.hla`
echo $hla_line
netMHCpan -BA -s -a $hla_line $indir/test/test.fasta/${name}.snp.fasta>$indir/test/netMHCpan/${name}.snp.out

less -S $indir/test/netMHCpan/${name}.snp.out|grep "Pos"|head -1 >$indir/test/netMHCpan/header
less -S $indir/test/netMHCpan/${name}.snp.out|grep "<= SB"|awk '{gsub("<= ","<=");print $0}'>$indir/test/netMHCpan/${name}.snp.pep
cat $indir/test/netMHCpan/header $indir/test/netMHCpan/${name}.snp.pep >$indir/test/netMHCpan/${name}.snp.sb
done

rm $indir/test/affinity.file
cd $indir/test/netMHCpan
ls $indir/test/netMHCpan/*.snp.sb|while read id 
do
checkfilesize()
{
          file=$(basename $id )
          name=${file%%.*}
          filename=${name}.snp.sb
          echo $name
          touch $indir/test/affinity.file
          filezise=`ls -l $filename|awk '{print $5}'`
          maxsize=$((0))  
          echo $filesize $maxsize                      #max memory
          if [ ! -s $filename ];then   
          echo "$name" >> $indir/test/affinity.file #backup
          rm  $filename
           fi
}
checkfilesize
done

echo -e "\n*** SKCM_netMHCpan_filter_indel.R started at $(date +'%T %F') ***\n"
Rscript $HOME/R/rscripts/SKCM_netMHCpan_filter.R
