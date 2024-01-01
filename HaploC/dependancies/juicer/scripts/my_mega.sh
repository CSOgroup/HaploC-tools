#!/bin/bash
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# MegaMap script. 
#
#                                                                       
# [topDir]    - Should contain the results of all base experiments
#
# From the top-level directory, the following two directories are created:
#                                                                              
# [topDir]/mega     - Location of result of processing the mega map
#
# Juicer version 1.5
# juicer_version="1.5.6" 

# usageHelp="Usage: ${0##*/} -g genomeID [-d topDir] [-s site] [-h]"
# genomeHelp="   genomeID must be defined in the script, e.g. \"hg19\" or \"mm10\" (default \"$genomeID\")"
# dirHelp="   [topDir] is the top level directory (default \"$topDir\") and must contain links to all merged_nodups files underneath it"
# siteHelp="   [site] must be defined in the script, e.g.  \"HindIII\" or \"MboI\" (default \"$site\"); alternatively, this can be the restriction site file"
# stageHelp="* [stage]: must be one of \"final\", \"postproc\", or \"early\".\n    -Use \"final\" when the reads have been combined into merged_nodups but the\n     final stats and hic files have not yet been created.\n    -Use \"postproc\" when the hic files have been created and only\n     postprocessing feature annotation remains to be completed.\n    -Use \"early\" for an early exit, before the final creation of the stats and\n     hic files"
# excludeHelp="   -x: exclude fragment-delimited maps from Hi-C mega map (will run much faster)"
# helpHelp="   -h: print this help and exit"

# printHelpAndExit() {
#     echo "$usageHelp"
#     echo "$genomeHelp"
#     echo "$dirHelp"
#     echo "$siteHelp"
#     echo "$stageHelp"
#     echo "$excludeHelp"
#     echo "$helpHelp"
#     exit "$1"
# }


while getopts "d:g:s:p:y:D:f:b:" opt; do
    case $opt in       
    g) genomeID=$OPTARG ;;
    d) topDir=$OPTARG ;;
    s) site=$OPTARG ;;
    p) genomePath=$OPTARG ;;  
    y) site_file=$OPTARG ;;
    D) juiceDir=$OPTARG ;;
    b) ligation=$OPTARG ;;
    [?]) printHelpAndExit 1;;
    esac
done



## Set ligation junction based on restriction enzyme
case $site in
    HindIII) ligation="AAGCTAGCTT";;
    DpnII) ligation="GATCGATC";;
    MboI) ligation="GATCGATC";;
    BglII) ligation="GATCGATC";;
    # BglII added by Yuanlong, 06-01-2019, reference: https://international.neb.com/tools-and-resources/selection-charts/recleavable-filled-5-overhangs
    none) ligation="XXXX";;
    *)  ligation="XXXX"
    site_file=$site
    echo "$site not listed as recognized enzyme, so trying it as site file."
    echo "Ligation junction is undefined";;
esac


echo $site_file
echo $ligation
echo $site


## Directories to be created and regex strings for listing files
megadir=${topDir}"/mega"
outputdir=${megadir}"/aligned"
tmpdir=${megadir}"/HIC_tmp"
export TMPDIR=${tmpdir}
outfile=${megadir}/lsf.out
#output messages
logdir="$megadir/debug"
logfile=${topDir}"/mega/log_mega.txt"

mkdir -p ${outputdir}

SECONDS=0
date > $logfile


mega_count=`find ${topDir}/mega/ -type f \( -name "merged_nodups.txt.gz" -o -name "merged_nodups.txt" \) | wc -l`
if [ "$mega_count" -ge "1" ]
then
    echo 'You already have merged_nodups.txt in mega dir. Remove to reprocess' > $logfile
    exit
fi

## Check for existing merged_nodups files:
merged_count=`find ${topDir} -type f \( -name "merged_nodups.txt.gz" -o -name "merged_nodups.txt" \) | wc -l`
if [ "$merged_count" -lt "1" ]
then
    echo "***! Failed to find at least one merged_nodups files under ${topDir}"
    exit 1
fi

echo n_merged_nodups file: $merged_count >> $logfile


## Yualong modified
merged_names=$(find ${topDir} -type f -name "merged_nodups.txt.gz" | awk '{print "<(gunzip -c",$1")"}' | tr '\n' ' ')
if [ ${#merged_names} -eq 0 ]
then
    merged_names=$(find ${topDir} -type f -name "merged_nodups.txt" | tr '\n' ' ')
fi
inter_names=$(find ${topDir} -type f -name "inter.txt" | tr '\n' ' ')

# echo $merged_names
# exit

# merged_names=$(find -L ${topDir} | grep merged_nodups.txt.gz | awk '{print "<(gunzip -c",$1")"}' | tr '\n' ' ')
# if [ ${#merged_names} -eq 0 ]
# then
#     merged_names=$(find -L ${topDir} | grep merged_nodups.txt | tr '\n' ' ')
# fi
# inter_names=$(find -L ${topDir} | grep inter.txt | tr '\n' ' ')



## Create temporary directory
if [ ! -d "$tmpdir" ]; then
    mkdir $tmpdir
    chmod 777 $tmpdir
fi

## Create output directory, used for reporting commands output
if [ ! -d "$logdir" ]; then
    mkdir "$logdir"
    chmod 777 "$logdir"
fi

## Arguments have been checked and directories created. Now begins
## the real work of the pipeline


    echo "I am here"
    echo ${merged_names} >> $logfile
    # Create top statistics file from all inter.txt files found under current dir
    awk -f ${juiceDir}/scripts/common/makemega_addstats.awk ${inter_names} > ${outputdir}/inter.txt

    sort -T ${tmpdir} -m -k2,2d -k6,6d ${merged_names} > ${outputdir}/merged_nodups.txt ## not merged
    echo "(-: Finished sorting all merged_nodups files into a single merge."
    
    # rm -r ${tmpdir}
    ${juiceDir}/scripts/common/statistics.pl -q 1 -o${outputdir}/inter.txt -s $site_file -l $ligation ${outputdir}/merged_nodups.txt

    # echo ${juiceDir}/scripts/juicer_tools pre -r 10000,20000,25000,40000,50000,100000,1000000 -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
    ${juiceDir}/scripts/common/juicer_tools pre -r 10000,20000,25000,40000,50000,100000,1000000 -s ${outputdir}/inter.txt -g ${outputdir}/inter_hists.m -q 1 ${outputdir}/merged_nodups.txt ${outputdir}/inter.hic ${genomeID}
 
    # gzip -f ${outputdir}/merged_nodups.txt

    echo " " >> $logfile
    echo "(-: Successfully completed making mega map. Done. :-)" >> $logfile

    duration=$SECONDS
    echo "$(($duration / 3600)) hours or $(($duration / 60)) minutes elapsed." >> $logfile
