#!/bin/bash

rm -r output
mkdir output
python cluster_assoc_p1.py $1 output/datalist1.txt
../../mscluster2/MSCluster_bin --list output/datalist1.txt --out-dir output/out1 --model-dir ../../mscluster2/Models --output-name out --assign-charges --sqs .1 --memory-gb 32
mkdir output/altered-mgf
python cluster_assoc_p2.py output/out1/mgf output/altered-mgf
python cluster_assoc_p1.py output/altered-mgf output/datalist2.txt
../../mscluster2/MSCluster_bin --list output/datalist2.txt --out-dir output/out2 --model-dir ../../mscluster2/Models --output-name out --assign-charges --sqs .1 --memory-gb 32
mkdir output/final-mgf
python cluster_assoc_p3.py output/out1/clust output/out2/clust output/out1/mgf output/datalist1.txt output/final-mgf
