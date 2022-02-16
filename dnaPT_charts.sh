#!/bin/bash
#
# Changelog
# V.0 | 02.11.22 - first version
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
   >>>      dnaPT_charts.sh           <<<
   **************************************
   Author: Clément Goubert - goubert.clement@gmail.com
   Last revision: 02/11/2022

   This script process an output folder of dnaPipeTE v1.3 and produces 3 graphs:
   - Barplot indicating the genomic proportion of each dnaPipeTE contig and the associated classification.
     A threshold -t (% genome) can be set such as contigs representing < t% of the genome will be groupped
     together.
   - Two piecharts: the first one represent the relative proportion of the different repeat categories, 
     and the second one indicate their proportion relative to the total genome.

   Note: to plot the landscape analysis, see the script "dnaPT_landscape.sh" included in this repositoty

   Dependencies:
   - R + package "ggplot2", "tidyverse" and "cowplot" (https://www.r-bloggers.com/2010/11/installing-r-packages/)

   ***************************************

   Usage: ./dnaPT_charts.sh -I <dataset_directory> [options]

   mendatory arguments:
    -I, --input-dir              dnaPipeTE output directory

   options:
    -p, --prefix                 prefix to append to the output filename: "<prefix>_charts.pdf"
    -t, --percent_threshold      barplot: repeats < -t % are groupped together as a single category (default 0.001%)
    -o, --output                 output folder (path); default: dnaPipeTE output directory
    -y, --ylim                   Max value for the y axis (genome %) [0-100]
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
   -p|--prefix)
       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
         PREF=$2
         shift 2
   #    else
   #      echo "Error: missing dataset B" >&2
   #      usage
   #      exit 1
       fi
       ;;    
    -y|--ylim)
       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
         YMAX=$2
         shift 2
   #    else
   #      echo "Error: missing prefix for dataset A" >&2
   #      usage
   #      exit 1
       fi
       ;;
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
   -t|--percent_threshold)
      if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
        TPERC=$2
        shift 2
      fi
      ;;
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
     # -U | --no-unknown)
     #     UNK=FALSE
     #     shift
     #    ;;
     # -S | --superfamily)
     #     SF=TRUE
     #     shift
     #    ;;         
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
PREF="${PREF:-dnaPipeTE}"
YMAX="${YMAX:-FALSE}"
TPERC="${TPERC:-0.001}"
# param check
echo "input dataset:      $DSA"
echo "output folder:      $OUTF"
echo "output prefix:      $PREF"
echo "custom ylim:        $YMAX"


############ START ############

mkdir -p $OUTF

Rscript - <<SCRIPT
################################################################################
# packages loading  / error if absent                                          #
################################################################################
packages <- c("ggplot2", "tidyverse", "cowplot")
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
rca<-read.table("$DSA/reads_per_component_and_annotation", fill = T, na.strings=c("","NA"))
rca<-separate(rca, V6, c("subclass", "superfamily"), sep = "/",fill = "right") 
rca\$superfamily[is.na(rca\$superfamily)]<-"Unknown"
rca\$subclass[is.na(rca\$subclass)]<-"Unknown"
# print("consolidate SF column...")
rca\$superfamily<-paste(rca\$subclass, rca\$superfamily, sep = "/")
# print("renames NA in SF columns...")
rca\$superfamily[rca\$superfamily == "Unknown/Unknown"]<-"Unknown"
names(rca)<-c("reads", "bp", "contig", "q_length", "target_name", "t_subc", "t_supf", "hit_cov")
reads.c<-as.numeric(system("grep -c '>' $DSA/renamed.blasting_reads.fasta", intern = T))
reads.b<-as.numeric(system("tail -n 2 $DSA/Counts.txt | head -n 1 | cut -f 2", intern = T))
threshold<-as.numeric("$TPERC")
rca<-rca[order(rca\$bp, decreasing = T),]
leftover<-c(
  sum(rca[rca\$bp < threshold*reads.b/100,]\$reads),
  sum(rca[rca\$bp < threshold*reads.b/100,]\$bp),
  rep(NA, 3),
  paste("repeats_under_", threshold, "%", sep = ""),
  paste("repeats_under_", threshold, "%", sep = ""),
  NA)
rca_t<-rca[rca\$bp >= threshold*reads.b/100,]
rca_t<-rbind(rca_t, leftover)
rca_t\$perc_g<-as.numeric(rca_t\$bp)/reads.b*100
rca_t\$x2<-cumsum(as.numeric(rca_t\$perc_g))
rca_t\$x1<-as.numeric(rca_t\$x2)-as.numeric(rca_t\$perc_g)

############ BARPLOT #############################
# pick the colors
cols<-read.table("$DIR/colors.land", sep = "\t")
cols[length(cols\$V1)+1,]<-c(paste("repeats_under_", threshold, "%", sep = ""), "grey10")
col.bars<-rep("", length((levels(as.factor(rca_t\$t_subc)))))
for(i in 1:length(levels(as.factor(rca_t\$t_subc)))){
  col.bars[i]<-cols\$V2[grep(pattern = paste("^", levels(as.factor(rca_t\$t_subc))[i], "$", sep = ""), x = cols\$V1)]
}

barplot<-ggplot(rca_t)+
  geom_rect(aes(xmin = x1, xmax = x2, ymin = 0, ymax = perc_g, fill = t_subc))+
  scale_fill_manual(values = col.bars, name = "repeat type")+
  xlab("cumulative genome %")+
  ylab("genome %")+
  theme(axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank())

############ PIE CHARTS #############################
rca_pie<-as.data.frame(tapply(rca_t\$perc_g, rca_t\$t_subc, sum))
rca_pie\$subc<-rownames(rca_pie)
names(rca_pie)<-c("perc_g", "subc")
rca_pie\$perc_r<-rca_pie\$perc_g/sum(rca_pie\$perc_g)*100


###### plot in % TE ######
df2 <- rca_pie %>% 
  mutate(csum = rev(cumsum(rev(perc_r))), 
         pos = perc_r/2 + lead(csum, 1),
         pos = if_else(is.na(pos), perc_r/2, pos))
df2\$subc<-factor(df2\$subc, levels = df2\$subc, labels = paste(paste0(round(df2\$perc_r,2), "% ",df2\$subc)))
pie_rel<-ggplot(df2, aes(x = "", y = perc_r, fill = subc))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y")+
  scale_fill_manual(values = col.bars, name = "repeat type (relative %)")+
  # geom_label_repel(data = df2,
  #                  aes(y = pos, label = paste0(round(perc_r,2), "%")),
  #                  size = 4.5, nudge_x = 1.5, show.legend = FALSE, fill = "white")+ #, col = col.bars)+
  theme_void()
  
###### plot in % genome #####
nr<-100-sum(rca_pie\$perc_g)
rca_pie_g<-rbind(rca_pie[,1:2], c(nr, "non-repetitive"))
col.bars_g<-c(col.bars, "grey90")
rca_pie_g\$perc_g<-as.numeric(rca_pie_g\$perc_g)
dfg <- rca_pie_g %>% 
  mutate(csum = rev(cumsum(rev(perc_g))), 
         pos = perc_g/2 + lead(csum, 1),
         pos = if_else(is.na(pos), perc_g/2, pos))
dfg\$subc<-factor(dfg\$subc, levels = dfg\$subc, labels = paste(paste0(round(dfg\$perc_g,2), "% ",dfg\$subc)))
pie_g<-ggplot(dfg, aes(x = "", y = perc_g, fill = subc))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y")+
  scale_fill_manual(values = col.bars_g, name = "repeat type (genome %)")+
  #geom_label_repel(data = dfg,
  #                 aes(y = pos, label = paste0(round(perc_g,2), "%")),
  #                 size = 4.5, nudge_x = 3,min.segment.length = 15, show.legend = FALSE, fill = "white")+ #, col = col.bars)+
  theme_void()

############ OUTPUT CHARTS ###########################
top<-plot_grid(pie_rel, pie_g, ncol = 2)
charts<-plot_grid(top, barplot, ncol = 1)

save_plot(
  file = paste("$PREF", "_charts.pdf", sep = ""),
  charts,
  path = "$OUTF",
  base_width = 10,
  base_height = 8
  )

SCRIPT