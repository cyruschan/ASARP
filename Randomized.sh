#!/bin/bash
rstart=21
rend=50
org_snvs=data/snp.list.fig3
rand_folder=random/test
rand_result=random/result
mkdir -p $rand_folder
mkdir -p $rand_result/tmp

echo "Generating random SNVs: $rstart to $rend"
echo perl -I random random/randAllSnvs.pl $org_snvs $rand_folder $rstart:$rend
perl -I random random/randAllSnvs.pl $org_snvs $rand_folder $rstart:$rend

