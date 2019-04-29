#!/bin/bash

FASTQDIR="/projects/ps-epigen/seqdata/"
export PATH=$PATH:/projects/ps-epigen/software/SNAP-CHIP_epicypher

## find snap log files
cd ~/data/outputs/libQCs/

find . -name "*snap*.log" | while read l
do
    pre=${l/.snap.log/};
    if [ ! -f ${pre}.snap.cnt.tab ]
    then
        fq=${FASTQDIR}/${pre##*/}
        echo ${pre}.snap.cnt.tab
        echo ${pre%/*}
        echo $fq
        bash epicypher.sh -i $fq -m 2 -k true -o ${pre%/*} 1 | tee ${pre}.snap.cnt.tab
    fi
done

## log all changes 
# find . -name "*snap*.cnt.tab" -mmin -30  > ~/data/logs/app/2019-04-08_snap_recover.txt
# cat  Set_245.txt |while read l;do  grep  ${l}/ ~/data/logs/app/2019-04-08_snap_recover.txt; done
