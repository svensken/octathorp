#! /bin/bash

octadir='/home/svensken/octathorp/'

# pass pdb numbers
pdb1='1WM1.pdb 142 242' # pose_from_pdb(pose, '1WM1.pdb')
                        # ERROR: too many tries in fill_missing_atoms!
                        # ERROR:: Exit from: src/core/conformation/Conformation.cc line: 2664
pdb2='2WUF.pdb 144 221'
pdb3='3ANS.pdb 366 471'
pdb4='3B12.pdb 133 222'
pdb5='3GZJ.pdb 119 197'


for host in "$pdb1" "$pdb2" "$pdb3" "$pdb4" "$pdb5" 
do
    set -- $host
    hostname=$1 # ****.pdb
    hres1=$2
    hres2=$3
    
    for guest in "$pdb1" "$pdb2" "$pdb3" "$pdb4" "$pdb5"
    do
        set -- $guest
        guestname=$1
        gres1=$2
        gres2=$3
        
        pair=$hostname'_'$guestname # ****.pdb_****.pdb
        mkdir -p $pair
        cd $pair

        waitforme=()
        for i in {1..30}
        do
            mkdir -p $i
            cd $i
            # /path/to/script.py /path/to/host.pdb hr1 hr2 /path/to/guest.pdb gr1 gr2 &
            # if pdb's pre-preped:
            # /path/to/script.py /path/to/host.pdb /path/to/guest.pdb &
            /home/svensken/octathorp/docky.py $octadir$host $octadir$guest &
            waitforme+=($!)
            cd -
        done
        
        wait ${waitforme[@]} # wait for this round of processes

        cd ..
    done

done


