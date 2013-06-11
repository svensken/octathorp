#!/bin/bash

ncpu=10

matchexe="/home/svensken/Rosetta/main/source/bin/match.linuxgccrelease"
statusupdate="/home/svensken/octathorp/status.update"


combopath="/home/svensken/octathorp/matching/combos"
for combodir in $combopath/*/
do
    cd $combodir
    id=${combodir: -10:4}
    flagdir="/home/svensken/octathorp/matching/flags/"$id
    flags="@$flagdir/general_match.flags @$flagdir/scaf_gfp.flags @$flagdir/substrate_gfp.flags"

    t1=$(date +%s)
    ls $combodir/*.pdb | xargs -n 1 -P $ncpu    $matchexe $flags 
    t2=$(date +%s)
    echo "$combodir $(expr $t2 - $t1)" >> $statusupdate
done





#for a in {1..1120}
#do
#    /home/svensken/Rosetta/main/source/bin/match.linuxgccrelease @../flags/general_match.flags @../flags/scaf_gfp.flags @../flags/substrate_gfp.flags -s /home/svensken/octathorp/patchdock/combos/3ANSbot.pdb_3ANScap.pdb.dir/output.txt.$a.pdb
#    echo $a > status.update
#done

#for b in {1..790}
#do
#    /home/svensken/Rosetta/main/source/bin/match.linuxgccrelease @../flags/general_match.flags @../flags/scaf_gfp.flags @../flags/substrate_gfp.flags -s /home/svensken/octathorp/patchdock/combos/3ANSbot.pdb_3B12cap.pdb.dir/output.txt.$b.pdb
#    total=$(expr 1120 + $b)
#    echo $total > status.update
#done




