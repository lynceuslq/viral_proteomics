#!/bin/bash

##########################here are the thing you need to define before running the ###################
export PATH="/PATH/TO/JAVA/bin/:$PATH" #do not use anything containing "TESTFAM", "TESTWORKDIR", "TESTINPATH" or "INPUTFASTA", the same in the lines below
INTERPROSCAN="/PATH/TO/INTERPROSCAN/interproscan-5.50-84.0/interproscan.sh"
fullnamelineage="/PATH/TO/LINEAGEFILE/fullnamelineage.dmp" # a file decompressed from new_dump.tar.gz obtained from NCBI, a faster way is extracting taxa derived from Orthornavirae and store it as Orthornavirae.fullnamelineage.dmp here


##########################you do not need to change anything below#####################
fam="TESTFAM"
wokingdir="TESTWORKDIR/$fam"
Infiledir="TESTINPATH"
inputf="INPUTFASTA"

echo -e "start to get seq for $fam at $(date)"

mkdir $wokingdir
grep -w "$fam" $fullnamelineage > $wokingdir/${fam// /}.fullnamelineage.dmp
cut -f1 $wokingdir/${fam// /}.fullnamelineage.dmp | while read acc; do awk -v a="$acc" '$3 == a' ${Infiledir// /}/${inputf// /}.allsig.list.alltaxa >> $wokingdir/${fam// /}.acclist; done

echo -e "starting seq extraction for $fam at $(date)"

cut -f1 $wokingdir/${fam// /}.acclist | sort | uniq | while read acc; do grep -w "$acc" -A1 ${Infiledir// /}/${inputf// /}.allsig.faa >> $wokingdir/${fam// /}.assigned.faa; done

echo -e "seq extraction completed for $fam at $(date)"

/hwfssz5/ST_INFECTION/GlobalDatabase/share/software/Miniconda3/envs.multi-user/cd-hit-4.8.1/bin/cd-hit -i $wokingdir/${fam// /}.assigned.faa -o $wokingdir/${fam// /}.clusetred.faa -c 0.9 -aS 0.8 -p 1

echo -e "seq clustering completed for $fam at $(date)"

echo -e "starting to scan $fam at $(date)"

$INTERPROSCAN -i $workingdir/${fam// /}.clusetred.faa -f tsv -dp -o $workingdir/${fam// /}.clusetred.intprsc

echo -e "finishing interproscan $fam at $(date)"

echo -e "start to make an archive of interproscan results on $fam at $(date)"
cat $workingdir/${fam// /}.clusetred.intprsc | cut -f1 | sort | uniq | while read acc;do echo -e "$(awk -v a="$acc" '$1 == a' $wokingdir/${fam// /}.acclist)\t|\t$fam\t|\t$(awk -v a="$acc" '$1 == a' $workingdir/${fam// /}.clusetred.intprsc | cut -f5 | tr "\n" ";")"; done > $workingdir/${fam// /}.clusetred.intprsc.archive

echo -e "start to generate a summary on interproscan results of $fam at $(date)"

ls $workingdir/${fam// /}.clusetred.intprsc | while read file; do echo $file; m=$( cut -f1 $file |sort | uniq | wc -l); echo $m; n=$( cut -f1 $file | cut -d "_" -f1,2 |sort | uniq | wc -l) ; echo $n; cut -f5 $file |sort | uniq | while read acc; do echo -e "$acc\t$(grep "$acc" $file | cut -f1 | sort | uniq | wc -l | awk -v s="$m" '{print $1 / s}')\t$(grep "$acc" $file | cut -f1 | sort | uniq | wc -l | awk -v s="$n" '{print $1 / s}')\t$(grep "$acc" $file | cut -f1 | sort | uniq | wc -l)\t$(grep "$acc" $file | cut -f4,6 | head -1)";done  | sort -n -r -k2,2 > $workingdir/${fam// /}.clusetred.intprsc.sum

ls $wokingdir/${fam// /}.clusetred.faa | cut -d "/" -f1 | sort | uniq | while read fam; do awk -v f="$fam" '{if($1 ~ /^>/) {print $0, "\t|\t", f} else {print $0} }' $wokingdir/${fam// /}.clusetred.faa > $wokingdir/${fam// /}.clustered.titlewithfam.faa ; done

echo -e "job completed with $fam at $(date)"
