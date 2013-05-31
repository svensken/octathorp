#! /bin/bash

octadir='/home/svensken/octathorp/'

# pass pose numbers
pdb1='1WM1.pdb' # 137 138'
pdb2='2WUF.pdb' # 173 192'
pdb3='3ANS.pdb' # 383 466'
pdb4='3B12.pdb' # 149 212'
pdb5='3GZJ.pdb' # 128 190'


for hostname in "$pdb1" #"$pdb1" "$pdb2" "$pdb3" "$pdb4" "$pdb5" 
do
    #set -- $host
    #hostname=$1 # ~~~~.pdb
    #hres1=$2
    #hres2=$3
    
    for guestname in "$pdb3" #"$pdb1" "$pdb2" "$pdb3" "$pdb4" "$pdb5"
    do
        #set -- $guest
        #guestname=$1
        #gres1=$2
        #gres2=$3
        
        pair=$hostname'_'$guestname # 1WM1.pdb_3GZJ.pdb
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
            /home/svensken/octathorp/docky.py $octadir/forkenny/1WM1bot_3ANScap_input.pdb & #$host $octadir$guest &
            waitforme+=($!)
            cd -
        done
        
        wait ${waitforme[@]} # wait for this round of processes

        cd ..
    done

done


