#!/bin/sh

hg19=../data/hg19.merged.to.ensg.all.tx.03.18.2011.txt

line=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
echo ">> Unzip $hg19 if necessary"
echo ""

echo "Checking existance of $hg19 ..."
if [ ! -e $hg19 ]; then
  echo "Unzipping $hg19.gz to $hg ..."
  gunzip -c $hg19.gz >$hg19
  if [ -e $hg19 ]; then
    echo "[OK]"
  else
    echo "[FAILED]"
    exit;
  fi
else 
  echo "[OK]"
fi

echo "See if the file size is reasonable"
sleep 2
ls -hs $hg19
