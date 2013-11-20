#!/bin/sh

echo ""
line=">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
demoresult="demo2.results"
echo $line
echo "ASARP demonstration at Starting Point 3: bedgraph and RNA SNV files are ready in demo2.data folder"
echo "Results will be put into folder $demoresult"
echo "NOTE: if you want the complete demo starting with SAM files, run demo.sh"
echo $line
rm -rf $demoresult
mkdir $demoresult
echo ""

echo ">> Check Statistics::R installation"
perl ../testR.pl
echo ""

echo ">> Start ASARP with demo2.config and demo.param (check docs for details)"
echo "perl -I ../ ../asarp.pl demoresults/asarp.results demo2.config demo.param >asarp.log"
echo ">>>> Check asarp2.log for the screen output when it is finished"
perl -I ../ ../asarp.pl $demoresult/asarp.results demo2.config demo.param >asarp.demo2.log
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
