#!/bin/sh

module load suhail python/3.8.5
cd outputs
rm *
../parse_mmcif_to_domains_new.py ../uploads/$1
for i in *
do
CODE=`echo $i | cut -d_ -f1`
FILENAME=${CODE}_summary.txt
pdbinfo $i | grep UNIP >> ${FILENAME}
done

echo summary file created

