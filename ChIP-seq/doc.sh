#!/bin/sh

# multiBigwigSummary bins -b ./coverage/CTCF_rep1.bw ./coverage/CTCF_rep2.bw ./coverage/NANOG_rep1.bw ./coverage/NANOG_rep2.bw \
# ./coverage/POU5F1_rep1.bw ./coverage/POU5F1_rep2.bw ./coverage/H3K27me3_rep1.bw ./coverage/H3K27me3_rep2.bw \
# ./coverage/H3K36me3_rep1.bw ./coverage/H3K36me3_rep2.bw ./coverage/H3K4me1_rep1.bw ./coverage/H3K4me1_rep2.bw \
# ./coverage/H3K4me3_rep1.bw ./coverage/H3K4me3_rep2.bw -o ChIP_all.npz
plotCorrelation -in ChIP_all.npz --corMethod pearson --whatToPlot heatmap -o ChIP_correlation.pdf --skipZeros \
--removeOutliers --outFileCorMatrix ChIP_all.tab --colorMap vlag --plotNumbers