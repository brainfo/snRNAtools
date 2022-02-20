#!/bin/sh

# cooler zoomify 4DNFIIZQDVOR.mcool::resolutions/10000 -p 16 -r 40000 --balance --balance-args --blacklist=/home/data4/jh/data/reference/Mus_musculus/encode/ENCFF547MET_flattened.bed --balance-args -f --balance-args --convergence-policy=error -o 4DNFIIZQDVOR_40k.mcool
# cooler zoomify 4DNFIWSPOUOW.mcool::resolutions/10000 -p 16 -r 40000 --balance --balance-args --blacklist=/home/data4/jh/data/reference/Mus_musculus/encode/ENCFF547MET_flattened.bed --balance-args -f --balance-args --convergence-policy=error -o 4DNFIWSPOUOW_40k.mcool
# cooler zoomify 4DNFIYEYHJUX.mcool::resolutions/10000 -p 16 -r 40000 --balance --balance-args --blacklist=/home/data4/jh/data/reference/Mus_musculus/encode/ENCFF547MET_flattened.bed --balance-args -f --balance-args --convergence-policy=error -o 4DNFIYEYHJUX_40k.mcool

cooltools diamond-insulation 4DNFIIZQDVOR_40k.mcool::resolutions/40000 1000000 > 4DNFIIZQDVOR_diamond_insulation.txt
cooltools diamond-insulation 4DNFIWSPOUOW_40k.mcool::resolutions/40000 1000000 > 4DNFIWSPOUOW_diamond_insulation.txt
cooltools diamond-insulation 4DNFIYEYHJUX_40k.mcool::resolutions/40000 1000000 > 4DNFIYEYHJUX_diamond_insulation.txt
cooltools diamond-insulation 4DNFIAUG25QU_40k.mcool::resolutions/40000 1000000 > 4DNFIAUG25QU_diamond_insulation.txt
