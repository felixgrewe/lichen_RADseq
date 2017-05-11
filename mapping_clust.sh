#!/bin/bash

mkdir -p clust.90
mkdir -p sam
FILES=/home/fgrewe/lichen_rad_corrected/3-pyrad_s2345/clust.90/*.consens.gz
for f in $FILES
do
  echo "Processing $f file..."
  name=${f%%.con*}
  name=${name##*/}
  #echo "Now $name for file..."
  bowtie2 -x /home/fgrewe/lichen_rad_corrected/reference/Rmela_ref -U $f --al-gz /home/fgrewe/lichen_rad_corrected/4-bowtie/clust.90/${name}aligned.consens.gz -S /home/fgrewe/lichen_rad_corrected/4-bowtie/sam/${name}.sam -f -N 1 --no-unal -L 20 -D 20 -R 3
  # take action on each file. $f store current file name
  # cat $f
done
