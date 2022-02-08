#!/bin/bash
#
# Changelog
# V.0 | 02.08.22 - first version
#
# Author: Clément Goubert - goubert.clement@gmail.com

###################################################################################
# PARSER: from https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f #
###################################################################################
# usage function
function usage()
{
   cat << HEREDOC

   **************************************
   >>>      dnaPT_landscapes.sh       <<<
   **************************************
   Author: Clément Goubert - goubert.clement@gmail.com
   Last revision: 02/08/2022

   This script perform at "TE landscape" analysis, i.e., it plots an histogram of the blastn divergence between raw reads (TE
   copies in the genomes) and their consensus sequences assembled in "Trinity.fasta". The script plots only putative TE seque
   -nce among "LINE", "SINE", "LTR", "DNA" and "Unknown" (a.k.a. "NA"). 

   Dependencies:
   - R + package "ggplot2" (https://www.r-bloggers.com/2010/11/installing-r-packages/)

   ***************************************

   Usage: ./dnaPT_landscape.sh -I <dataset_directory> [options]

   mendatory arguments:
    -I, --input-dir              dnaPipeTE output directory

   options:
    -o, --output                 output folder (path); default: dnaPipeTE output directory
    -S, --subclass               Plot with subclass information (instead of Class)
    -h, --help                   Prints this message and exit

HEREDOC
} 

# if no parameter given, output help and qui
if [[ $# -eq 0 ]] ; then
    echo '   **********************************'
    echo '   Error! No mendatory argument given'
    echo '   **********************************'
    usage
    exit 0
fi

# parameters parser
PARAMS=""
while (( "$#" )); do
  case "$1" in
# flags with arguments
    -I|--input-dir )
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        DSA=$2
        shift 2
      else
        echo "Error: missing dnaPipeTE directory" >&2
        usage
        exit 1
      fi
      ;;
   # -B|--dir_B)
   #    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
   #      DSB=$2
   #      shift 2
   #    else
   #      echo "Error: missing dataset B" >&2
   #      usage
   #      exit 1
   #    fi
   #    ;;    
   # -a|--pref_A)
   #    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
   #      PREFA=$2
   #      shift 2
   #    else
   #      echo "Error: missing prefix for dataset A" >&2
   #      usage
   #      exit 1
   #    fi
   #    ;;
   # -b|--pref_B)
   #    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
   #      PREFB=$2
   #      shift 2
   #    else
   #      echo "Error: missing prefix for dataset B" >&2
   #      usage
   #      exit 1
   #    fi
   #    ;;       
   -o|--output)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        OUTF=$2
        shift 2
        else
        echo "Error: missing output folder" >&2
        usage
        exit 1
      fi
      ;;
    # -p|--percent_threshold)
    #    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
    #      TPERC=$2
    #      shift 2
    #    fi
    #    ;;
    #  -e|--ecp_threshold)
    #    if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
    #      TECP=$2
    #      shift 2
    #    fi
    #    ;;
#     -F | --full-length-alpha)
#       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#         FULL_ALPHA=$2
#         shift 2
#       fi
#       ;;   
#     -y | --auto-y)
#       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#         AUTO_Y=$2
#         shift 2
#       fi
#       ;;  
#     -m | --min-orf)
#      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#         MINORF=$2
#         shift 2
#       fi
#       ;; 
     -h | --help)
        usage
        exit 1
        ;;
#    -o | --output)
#      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#        OUTPUT=$2
#        shift 2  
#      fi
#       ;;
# boolean flags            
     #  -T | --te_only)
     #      TE=TRUE
     #      shift
     #    ;; 
     # -E | --ecp)
     #     ECP=TRUE
     #     shift
     #    ;;
     -S | --subclass)
         SUB=TRUE
         shift
        ;;         
    -*|--*=) # unsupported flags
      echo "Error: Unsupported argument $1" >&2
      usage
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac # <- end of case
done
# set positional arguments in their proper place
eval set -- "$PARAMS"

####################################################################
# MAIN:                                                            #
####################################################################

# get script launch dir, from https://stackoverflow.com/a/246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
# asign default value and print parameters
OUTF="${OUTF:-$DSA}"
SUB="${SUB:-FALSE}"

# param check
echo "input dataset:      $DSA"
echo "output folder:      $OUTF"
echo "Subclass level:     $SUB"

############ START ############

mkdir -p $OUTF

Rscript - <<SCRIPT
################################################################################
# packages loading  / error if absent                                          #
################################################################################
packages <- c("ggplot2")
# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  print(paste("ERROR: the following R packages are missing: ", packages[!installed_packages], sep = ""))
  print("quitting...")
  quit(status=1)
  #install.packages(packages[!installed_packages])
}
# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

################################################################################
# MAIN                                                                         #
################################################################################
# read main table and shell variables
print("load variables...")
# join the reads with annotations and format table for R
land<-read.table(sep = "\t", text=system("join -a1 -12 -21 -o 1.3,2.4,2.5  $DSA/Annotation/sorted_blast3 $DSA/Annotation/one_RM_hit_per_Trinity_contigs | awk '/LINE/ { print \$0 \"\\t\" \$3; next} /LTR/ {print \$0 \"\\t\" \$3; next} /SINE/ {print \$0 \"\\tSINE\"; next} /DNA/ {print \$0 \"\\tDNA\"; next} /MITE/ {print \$0 \"\\tMITE\";next} {if (NF == 3) {print \$0\"\tOthers\"} else {print \$0\"\\tNA\\tUnknown\\tUnknown\"}}' | sed 's/ /\t/g;s/\t\t\t/\\t/g' ", intern = T))
reads.c<-as.numeric(system("grep -c '>' $DSA/renamed.blasting_reads.fasta", intern = T))

# pick the colors
cols<-read.table("$DIR/colors.land", sep = "\t")
col.lands<-rep("", length((levels(as.factor(land\$V4)))))
for(i in 1:length(levels(as.factor(land\$V4)))){
  col.lands[i]<-cols\$V2[grep(pattern = paste("^", levels(as.factor(land\$V4))[i], "$", sep = ""), x = cols\$V1)]
}

# plot
lscapes<-ggplot(land, aes(100-V1, fill = V4))+
  geom_histogram(aes(y=..count../reads.c*100), binwidth = 1)+
  scale_fill_manual(values = col.lands)+
  ylab("genome %")+
  xlab("blastn divergence (read vs dnaPipeTE contig)")

# export
ggsave(
    file = "landscapes.pdf",
    plot = lscapes,
    device = "pdf",
    path = "$OUTF",
    scale = 1,
    width = 2600,
    height = 2000,
    units = "px"
   )

SCRIPT

echo "All done, results in $OUTF/"