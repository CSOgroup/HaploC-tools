#!/usr/bin/python

# author: Peter Edge + Y.LIU
# 12/19/2016
# ref: https://github.com/vibansal/HapCUT2/blob/master/recipes/HiC_10X/Snakefile

import os,sys

input_mate1 = sys.argv[1]
input_mate2 = sys.argv[2]
output_repaired = sys.argv[3]
HAPCUT2_REPO = sys.argv[4]

################################################################################

utilities_dir = os.path.join(HAPCUT2_REPO, 'utilities')
sys.path.append(utilities_dir)
import HiC_repair

################################################################################

HiC_repair.repair_chimeras(input_mate1, input_mate2, output_repaired, min_mapq=10)
