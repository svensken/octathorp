#!/bin/bash

ncpu=10

matchexe="/home/svensken/Rosetta/main/source/bin/match.linuxgccrelease"
statusupdate="/home/svensken/octathorp/status.update"


combopath="/home/svensken/octathorp/matching/combos"
for combodir in $combopath/3ANS*/ #$combopath/*/
do
    cd $combodir
    id=${combodir: -10:4}
    flagdir="/home/svensken/octathorp/matching/flags/"$id
    flags="@$flagdir/general_match.flags @$flagdir/scaf_gfp.flags @$flagdir/substrate_gfp.flags"

    t1=$(date +%s)
    ls $combodir/*.pdb | xargs -n 1 -P $ncpu    $matchexe $flags -s #arg
    t2=$(date +%s)
    echo "$combodir $(expr $t2 - $t1)" >> $statusupdate
done




