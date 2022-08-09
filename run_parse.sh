#!/bin/sh

# Run AFsegment
module load suhail python/3.8.5
rm public/*.jpg
cd outputs
rm *
../parse_mmcif_to_domains_new.py ../uploads/$1

# Create summary of domain info
for i in *
do
CODE=`echo $i | cut -d_ -f1`
FILENAME=${CODE}_summary.txt
pdbinfo $i | grep UNIP >> ${FILENAME}
done
echo summary file created

# Python process domain info and create .pml file
FILENAME=`find *summary.txt`
python ../summary_to_pml.py $FILENAME

# Create PyMol Image
module load pymol
pymol ../uploads/$1 -d zoom -c colbydom.pml -g colbydom.png
convert -quality 75% colbydom.png ../public/colbydom.jpg

rm ../uploads/*

# get ../af_files_cif/AF-P16114-F1-model_v2.cif example_cif/.