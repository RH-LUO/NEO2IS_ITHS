#!/bin/bash
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -p cn
source $HOME/miniconda3/bin/activate
indir=$HOME/lrh/SKCM/Neoantigen/SKCM.indel_results
outdir=$HOME/lrh/SKCM/Neoantigen/SKCM.indel_results/fastaFiles

rm $indir/fasta.missed.file
cd $outdir
ls $outdir/*.reformat.fasta|while read id
do
file=$(basename $id )
tumor=${file%%.*}
echo $file $tumor

checkfilesize()
{

       time=$(date "+%Y-%m-%d %H:%M:%S")   #
       filename=${tumor}.reformat.fasta
       echo $id $filename
       #if [ ! -f $filename ];then          #
             # touch $outdir/$filename
              touch $indir/fasta.missed.file
             # echo "" >> $outdir/$filename
          #    echo "" >> ../fasta.log
        #  fi
          filezise=`ls -l $filename|awk '{print $5}'`
          maxsize=$((0))  
          echo $filesize $maxsize                      #
       if [ ! -s $filename ];then   
         echo "$tumor" >> $indir/fasta.missed.file #
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

outdir=$HOME/lrh/SKCM/Neoantigen/SKCM.indel_results/test/test.fasta
cd $outdir

echo -e "\n*** SKCM_neopeptide_filter_indel.R started at $(date +'%T %F') ***\n"
Rscript $HOME/R/rscripts/SKCM_neopeptide_filter_indel.R

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
          echo $filesize $maxsize                      #
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
cat $indir/test/test.fasta/${name}.test.fasta|grep -v "silent"|grep -v "to )"|grep -v "to NA)"|awk '{gsub(")",")\n"); print $0}'>$indir/test/test.fasta/${name}.indel.fa
 done

rm $indir/test/nonefa.file
cd $indir/test/test.fasta
ls $indir/test/test.fasta/*.indel.fa |while read id 
do
checkfilesize()
{
          file=$(basename $id )
          name=${file%%.*}
          filename=${name}.indel.fa
          echo $name
          touch $indir/test/nonefa.file
          filezise=`ls -l $filename|awk '{print $5}'`
          maxsize=$((0))  
          echo $filesize $maxsize                      #
          if [ ! -s $filename ];then   
          echo "$name" >> $indir/test/nonefa.file #
          rm  $filename
  #        rm $outdir/pos/${name}.indel.bed
           fi
}
checkfilesize
done

ls $indir/test/test.fasta/*.indel.fa|while read id
do
file=$(basename $id )
name=${file%%.*}
echo $name  
cat  $indir/test/test.fasta/${name}.indel.fa|grep ">" >$indir/test/pos/${name}.indel.sites

cat $indir/test/pos/${name}.indel.sites| awk  'BEGIN { FS = "position" } ; { print $2 }'|awk  'BEGIN { FS = "-" }; { print $1 }'|awk '{print $1}' >$indir/test/pos/${name}.indel.pos

paste  $indir/test/pos/${name}.indel.pos $indir/test/pos/${name}.indel.sites > $indir/test/pos/${name}.indel.position

cat $indir/test/pos/${name}.indel.position|awk 'BEGIN{OFS=FS=" "}{if($8!="has") {tss= $1; end=$12; tss_up=tss"+end"; tss_dw=tss-10;} else {tss= $1; end=$10; tss_up=tss+10; tss_dw=tss-10;} if(tss_dw<0) tss_dw=0;print $2, tss, tss_dw, tss_up, end;}'|awk '{gsub(">","");print $0}'|awk '{gsub(")","");print $0}'>$indir/test/pos/${name}.indel.bed
done

echo -e "\n*** SKCM_neopeptide_shorten_indel.R started at $(date +'%T %F') ***\n"
Rscript $HOME/R/rscripts/SKCM_neopeptide_shorten_indel.R

ls $indir/test/test.fasta/*.indel.fa|while read id
do
file=$(basename $id )
name=${file%%.*}
echo $name
hla_line=`cat $HOME/lrh/SKCM/Neoantigen/SKCM.snp_results/test/hla/${name}.hla`
if [ ! -e $indir/test/netMHCpan/${name}.indel.sb ]; then  
echo $hla_line

netMHCpan -BA -s -a $hla_line $indir/test/test.fasta/${name}.indel.fasta>$indir/test/netMHCpan/${name}.indel.out

less -S $indir/test/netMHCpan/${name}.indel.out|grep "Pos"|head -1 >$indir/test/netMHCpan/header
less -S $indir/test/netMHCpan/${name}.indel.out|grep "<= SB"|awk '{gsub("<= ","<=");print $0}'>$indir/test/netMHCpan/${name}.indel.pep
cat $indir/test/netMHCpan/header $indir/test/netMHCpan/${name}.indel.pep >$indir/test/netMHCpan/${name}.indel.sb
fi
done

rm $indir/test/affinity.file
cd $indir/test/netMHCpan
ls $indir/test/netMHCpan/*.indel.sb|while read id 
do
checkfilesize()
{
          file=$(basename $id )
          name=${file%%.*}
          filename=${name}.indel.sb
          echo $name
          touch $indir/test/affinity.file
          filezise=`ls -l $filename|awk '{print $5}'`
          maxsize=$((0))  
          echo $filesize $maxsize                      #
          if [ ! -s $filename ];then   
          echo "$name" >> $indir/test/affinity.file #
          rm  $filename
           fi
}
checkfilesize
done

echo -e "\n*** SKCM_netMHCpan_filter_indel.R started at $(date +'%T %F') ***\n"
Rscript $HOME/R/rscripts/SKCM_netMHCpan_filter_indel.R
