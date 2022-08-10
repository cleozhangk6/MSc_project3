#!/bin/sh

# Run AFsegment
module load python/3.8.5
cd public/outputs
rm *
../segment_alphafold.py ../uploads/$1

# Create summary of domain info
rm temp.pdb
# for i in *
# do
# CODE=`echo $i | cut -d_ -f1`
# FILENAME=${CODE}_summary.txt
# pdbinfo $i | grep UNIP >> ${FILENAME}
# done
# echo summary file created
for i in *
do
CODE=`echo $i | cut -d_ -f1`
FILENAME=${CODE}_summary.txt
echo $i | cut -d. -f1 >> ${FILENAME}
done
echo summary file created without pdbinfo

# Python process domain info and create .pml file
FILENAME=`find *summary.txt`
python ../summary_to_pml.py $FILENAME $1
echo jmol spt file $FILENAME $1 created

# # Create PyMol Image
# module load pymol
# pymol ../uploads/$1 -d zoom -c colbydom.pml -g colbydom.png
# convert -quality 75% colbydom.png ../public/colbydom.jpg
