#!/bin/sh

line=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

result="preproc"
rm -rf $result # clear any old result

mkdir -p $result/R1
mkdir -p $result/R2
mkdir -p $result/merged

echo $line
echo "ASARP demonstration pipeline: Preprocessing Part"
echo "Results will be put into folder $result"
echo $line
echo ""

echo ">> Check Statistics::R installation"
perl ../testR.pl
echo ""

echo ">> ASARP Preprocessing"
echo ">>>> demo SAM files from 2 replicates (R1, R2) in data folder:"
ls data/chr*.sam
echo ""

echo ">>>> Remove PCR duplicates: rmDup.pl (or samtools):"
for rep in 1 2
do
  echo ">>>> R$rep (check rmDup.R$rep.log for the screen output)"
  rlogfile=rmDup.R$rep.log
  rm -rf $rlogfile # clear the old log file
  #echo ">>>> USAGE: perl ../rmDup.pl input_sam_file output_sam_file is_paired_end"
  for chr in 1 5 10
  do
    #you can print out the command by uncommenting the line below
    #echo "perl -I ../ ../rmDup.pl data/chr$chr.R$rep.sam $result/chr$chr.R$rep.rmDup.sam 1 >>$rlogfile"
    perl -I ../ ../rmDup.pl data/chr$chr.R$rep.sam $result/R$rep/chr$chr.rmDup.sam 1 >>$rlogfile
    echo "chr$chr DONE" >> $rlogfile
  done
done
echo ">>>> Merge the SAM files: mergeSam.pl (or samtools):"
echo ">>>> A replicate folder list is needed, as shown in rep.demo.lst:"
cat rep.demo.lst
mlogfile=mergeSam.log
echo ">>>> results will be in $result/merged folder"
#echo ">>>> Start merging (check $mlogfile for the screen output)"
#rm -rf $mlogfile # clear the old log file 
echo ">>>> merging chr*.rmDup.sam in $result/R1 $result/R2 (check $mlogfile for the screen output)"
#you can print out the command by uncommenting the line below
#echo "perl -I ../ ../mergeSam.pl rep.demo.lst chr rmDup.sam $result/merged >$mlogfile"
perl -I ../ ../mergeSam.pl rep.demo.lst chr rmDup.sam $result/merged >$mlogfile
echo ""

if [ -e $result/merged/chr10.sam ]; then
  echo ">> Duplicate removal and merging FINISHED."
  echo ""
else
  echo ""
  echo "!! $result/merged/chr10.sam not generated. Duplicate removal and merging FAILED. Check step logs and error messages."
  exit;
fi

echo ">>>> Process reads to get SNV allelic reads and bedgraph files: procReads.pl"
echo "     (or samtools for SNV calling, and bedtools coverage -split for bedgraph generation)"
snvlist="dna.snv.demo.lst"
echo ">>>> A heterozyous snv list is needed, either from dbSNP or genomic sequencing. Example: $snvlist (from genomic sequencing)"
head -2 $snvlist
echo "..."
echo "The demo data are paired-end strand-specific where pair 2 is sense, so the setting is: is_paired_end=1 and is_strand_sp=2"
echo ">>>> Start processing... check procReads.chr*.log for the screen output"
for chr in 1 5 10
do
  #echo "check procReads.chr$chr.log for the screen output"
  #echo "perl -I ../ ../procReads.pl chr$chr $result/merged/chr$chr.sam $snvlist $result/chr$chr.snv $result/chr$chr.bed 1 2 \"demo track chr$chr\" >procReads.chr.$chr.log"
  perl -I ../ ../procReads.pl chr$chr $result/merged/chr$chr.sam $snvlist $result/chr$chr.snv $result/chr$chr.bed 1 2 "demo chr$chr" >procReads.chr$chr.log
done
echo ""

rnasnvlist="rna.snv.demo.lst"
echo ">>>> Merge all SNVs with allelic reads (a.k.a RNA SNVs) and generate plots"
echo "Because the SNVs are from genomic sequencing, mono=0 is used to keep all mono-allelic SNVs"
echo "For SNVs called using RNA-Seq data only, one can filter mono-alleleic SNVs < x via setting mono=x"
echo "Results will be output to $result/$rnasnvlist (check mergeSnvs.log for the screen output)"
perl -I ../ ../mergeSnvs.pl $result/ snv mono=0 $result/$rnasnvlist 1 >mergeSnvs.log
echo ""

if [ -e $result/$rnasnvlist ]; then
  echo ">> $result/$rnasnvlist generated. Preprocessing Part ALL FINISHED."
else
  echo "!! $result/$rnasnvlist not generated. Preprocessing Part FAILED. Check step logs and error messages."
  exit;
fi

echo ""
demoresult="demo.results"
echo $line
echo "ASARP demonstration pipeline: ASARP Main Part"
echo "Results will be put into folder $demoresult"
echo $line
rm -rf $demoresult
mkdir $demoresult
echo ""

echo ">> Start ASARP with demo.config and demo.param (check docs for details)"
echo "perl -I ../ ../asarp.pl demoresults/asarp.results demo.config demo.param >asarp.log"
echo ">>>> Check asarp.log for the screen output when it is finished"
perl -I ../ ../asarp.pl $demoresult/asarp.results demo.config demo.param >asarp.log
echo ""

if [ -e $demoresult/asarp.results ]; then
  echo ">> ASARP FINISHED. Summary from $demoresult/asarp.results"
  head -13 $demoresult/asarp.results
  echo ">>>> check $demoresult/asarp.results.*.prediction for more detailed predictions"
  echo "Demo SUCCESSFUL"
else
  echo "!! $demoresult/asarp.results not generated. ASARP FAILED. Check step logs and error messages."
  exit;
fi


