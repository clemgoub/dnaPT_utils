#!/bin/bash
#
# Changelog
# V.0 | 02.04.22 - first version
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

   This script search and measure the relative abundance of shared TE families between two datasets analyzed with dnaPipeTE             
   
   Dependencies:
   - CD-HIT
   - R                                      
   ***************************************

   Usage: ./dnaPT_compare.sh -A <dataset_A_directory> -a <prefix_A> -B <dataset_B_directory> -b <prefix_B> -o <output_folder>

   mendatory arguments:
    
    -A, --dir_A                  dnaPipeTE output directory for dataset A (path)
    -B, --dir_B                  dnaPipeTE output directory for dataset B (path)
    -a, --pref_a                 prefix for dataset A (string)
    -b, --pref_b                 prefix for dataset B (string)
    -o, --output                 output folder (path)

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
    -A|--dir_A)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        DSA=$2
        shift 2
      else
        echo "Error: missing dataset A" >&2
        usage
        exit 1
      fi
      ;;
   -B|--dir_B)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        DSB=$2
        shift 2
      else
        echo "Error: missing dataset B" >&2
        usage
        exit 1
      fi
      ;;    
   -a|--pref_A)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        PREFA=$2
        shift 2
      else
        echo "Error: missing prefix for dataset A" >&2
        usage
        exit 1
      fi
      ;;
   -b|--pref_B)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        PREFB=$2
        shift 2
      else
        echo "Error: missing prefix for dataset B" >&2
        usage
        exit 1
      fi
      ;;       
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
#    -f|--full-length-threshold)
#       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#         FL=$2
#         shift 2
#       fi
#       ;;
#     -a|--alpha)
#       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#         ALPHA=$2
#         shift 2
#       fi
#       ;;
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
#     -h | --help)
#        usage
#        exit 1
#        ;;
#    -o | --output)
#      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
#        OUTPUT=$2
#        shift 2  
#      fi
#       ;;
# # boolean flags            
#     -t | --tables)
#         TABLES=TRUE
#         shift
#       ;; 
#     -T | --all-Tables)
#         ALLTAB=TRUE
#         shift
#        ;;
#     -D | --emboss-dotmatcher)
#         DOTMA=TRUE
#         shift
#        ;;         
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
OUTF="${OUTF:-compa}"

# param check
echo "dataset A:     $DSA"
echo "dataset B:     $DSB"
echo "prefix A:      $PREFA"
echo "prefix B:      $PREFB"
echo "output folder: $OUTF"

mkdir -p $OUTF

# combine the contigs from each species’ dnaPipeTE run + add species name
#fileA=$
#fileB=$
cat <(sed -E 's/>/>'"$PREFA"_'/g' "$DSA"/Trinity.fasta) <(sed -E 's/>/>'"$PREFB"_'/g' "$DSB"/Trinity.fasta) > $PREFA''_$PREFB''_dnaPipeTE_contigs.fasta

# cluster sequences using CD-HIT-EST
cd-hit-est -i $PREFA''_$PREFB''_dnaPipeTE_contigs.fasta -o $PREFA''_$PREFB -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -bak 1 -T 8

# clean up cd-hit outputs before joining conts for each species
sort -k1,1n $PREFA''_$PREFB.bak.clstr | sed 's/>//g;s/nt,//g;s/\.\.\.//g;s/\*/REP/g;s/at//g' | awk '/REP/ {print $1"\t"$2"\t"$3"\t"$4} !/REP/ {print $1"\t"$2"\t"$3"\tin_cluster"}'> $PREFA''_$PREFB''.clean.clstr

# gather total bp sampled per species
AC=$(grep 'Total' "$DSA"/Counts.txt | tail -n 1 | cut -f 2)
BC=$(grep 'Total' "$DSB"/Counts.txt | tail -n 1 | cut -f 2)

# joint annotations and counts
join -a1 -13 -21 <(sort -k3,3 $PREFA''_$PREFB''.clean.clstr) <(cat <(awk -v count="$AC" -v prefA="$PREFA" '{print prefA"_"$3"\t"$5"\t"$6"\t"$1"\t"$2"\t"count}' $DSA/reads_per_component_and_annotation) <(awk -v count="$BC" -v prefB="$PREFB" '{print prefB"_"$3"\t"$5"\t"$6"\t"$1"\t"$2"\t"count}' $DSB/reads_per_component_and_annotation) | sort -k1,1) | awk '{if (NF == 7) {print $1"\t"$2"\t"$3"\t"$4"\tNA\tNA\t"$5"\t"$6"\t"$7} else if (NF == 4) {print $1"\t"$2"\t"$3"\t"$4"\tNA\tNA\tNA\t0\t1"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' | sort -k 2,2n -k1,1  > $PREFA''_$PREFB''_R.tsv
