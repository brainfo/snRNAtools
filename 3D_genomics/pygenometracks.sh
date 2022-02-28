#!/bin/sh

#cat Agglo_Cluster1.bed9 Agglo_Cluster2.bed9 Agglo_Cluster3.bed9 | sort -V -k1,1 -k2,2 - > all_tad.bed9

pyGenomeTracks --tracks explain.ini --region chr1:16000000-26000000 -o chr1_16M_26M.pdf
