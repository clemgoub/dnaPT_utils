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
    -p, --percent_threshold      min. percent genome to plot (default = 0) / not used if -E/-ecp selected
    -E, --ecp                    perform comparison in equivalent copy (dataset counts in bp / representative sequence size)
    -e, --ecp_threshold          min. equivalent copy/ies to plot (default = 0)

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
    -p|--percent_threshold)
       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
         TPERC=$2
         shift 2
       fi
       ;;
     -e|--ecp_threshold)
       if [ -n "$2" ] && [ ${2:0:1} != "-" ]; then
         TECP=$2
         shift 2
       fi
       ;;
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
# boolean flags            
     # -P | --percent)
     #     PERC=TRUE
     #     shift
     #   ;; 
     -E | --ecp)
         ECP=TRUE
         shift
        ;;
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
#PERC="${PERC:-TRUE}"
ECP="${ECP:-FALSE}"
TPERC="${TPERC:-0}"
TECP="${TECP:-0}"

# param check
echo "dataset A:          $DSA"
echo "dataset B:          $DSB"
echo "prefix A:           $PREFA"
echo "prefix B:           $PREFB"
echo "output folder:      $OUTF"
#echo "percent:            $PERC"
echo "ecp:                $ECP"
echo "percent threshold:  $TPERC"
echo "ecp threshold:      $TECP"

# check exclusive parameters
# if [ ${PERC} == TRUE ] && [ ${ECP} == TRUE ]; then
#   echo "ERROR: requires only one of -P (percent) or -E (equivalent copy)"
#   usage
#   exit 1
# fi

############ START ############

mkdir -p $OUTF

# check if clustering done
if [ -s $OUTF/$PREFA''_$PREFB''.clean.clstr ]
  then
   echo "clustering files found and not empty, skipping to R analysis..."
  else
   cat <(sed -E 's/>/>'"$PREFA"_'/g' "$DSA"/Trinity.fasta) <(sed -E 's/>/>'"$PREFB"_'/g' "$DSB"/Trinity.fasta) > $OUTF/$PREFA''_$PREFB''_dnaPipeTE_contigs.fasta

   # cluster sequences using CD-HIT-EST
   cd-hit-est -i $OUTF/$PREFA''_$PREFB''_dnaPipeTE_contigs.fasta -o $OUTF/$PREFA''_$PREFB -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -bak 1 -T 8

   # clean up cd-hit outputs before joining conts for each species
   sort -k1,1n $OUTF/$PREFA''_$PREFB.bak.clstr | sed 's/>//g;s/nt,//g;s/\.\.\.//g;s/\*/REP/g;s/at//g' | awk '/REP/ {print $1"\t"$2"\t"$3"\t"$4} !/REP/ {print $1"\t"$2"\t"$3"\tin_cluster"}'> $OUTF/$PREFA''_$PREFB''.clean.clstr

   # gather total bp sampled per species
   AC=$(grep 'Total' "$DSA"/Counts.txt | tail -n 1 | cut -f 2)
   BC=$(grep 'Total' "$DSB"/Counts.txt | tail -n 1 | cut -f 2)

   # joint annotations and counts
   join -a1 -13 -21 <(sort -k3,3 $OUTF/$PREFA''_$PREFB''.clean.clstr) <(cat <(awk -v count="$AC" -v prefA="$PREFA" '{print prefA"_"$3"\t"$5"\t"$6"\t"$1"\t"$2"\t"count}' $DSA/reads_per_component_and_annotation) <(awk -v count="$BC" -v prefB="$PREFB" '{print prefB"_"$3"\t"$5"\t"$6"\t"$1"\t"$2"\t"count}' $DSB/reads_per_component_and_annotation) | sort -k1,1) | awk '{if (NF == 7) {print $1"\t"$2"\t"$3"\t"$4"\tNA\tNA\t"$5"\t"$6"\t"$7} else if (NF == 4) {print $1"\t"$2"\t"$3"\t"$4"\tNA\tNA\tNA\t0\t1"} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' | sort -k 2,2n -k1,1  > $OUTF/$PREFA''_$PREFB''_R.tsv
fi

Rscript - <<SCRIPT
################################################################################
# packages loading  / error if absent                                          #
################################################################################
packages <- c("ggplot2", "gridExtra", "tidyr", "reshape2")
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
input<-paste("$OUTF", "/", "$PREFA", "_", "$PREFB", "_R.tsv", sep = "") # maybe change such as the file is a variable already?
data1<-"$PREFA"
data2<-"$PREFB"
out<-"$OUTF"
ecp_status<-"$ECP"
pc_T<-as.numeric("$TPERC")
ecp_T<-as.numeric("$TECP")

print("load table...")
DD<-read.table(input)
# create a column with dataset name
print("split column...")
DD<-separate(DD, V1, c("V1.1", "V1.2"), sep = "_comp_") 
# rename the columns
print("rename columns...")
names(DD)<-c("dataset", "contig", "cluster", "length", "status", "TE_name", "TE_Class", "reads", "bp", "total_SP")
# create a column with the percentage genome per contig per dataset
print("compute percentages...")
DD\$pc<-DD\$bp/DD\$total*100 
# aggregate the count (in %) for each cluster and each dataset
print("run dcast and convert table...")
counts<-as.data.frame(t(dcast(DD, formula = dataset~as.factor(cluster), value.var = "pc", fun.aggregate = sum))) 
counts\$cluster<-rownames(counts)
counts<-counts[-1,]
names(counts)<-c(data1, data2, "cluster")
# create a table with only one annotation per cluster (from REPresentative sequence of CD-HIT)
print("filter down to shared clusters...")
Drep<-DD[DD\$status == "REP",]
counts\$TE_class<-Drep\$TE_Class
#print(head(counts))

print("new additions...")
print("dcast bp...")
bp<-as.data.frame(t(dcast(DD, formula = dataset~as.factor(cluster), value.var = "bp", fun.aggregate = sum)))
bp<-bp[-1,]
print("rename column with paste function...")
names(bp)<-c(paste(data1, "_bp", sep = ""), paste(data2, "_bp", sep = ""))
print("cbind tables...")
counts<-cbind(counts, bp)
counts\$REPlen<-Drep\$length
print("get rep name...")
counts\$REPname<-paste(Drep\$dataset, Drep\$contig, sep = "_")
print("split class...")
counts<-separate(counts, TE_class, c("Class", "Super_family"), sep = "/",fill = "right") 
print("rename Unknowns...")
counts\$Class[is.na(counts\$Class)]<-"Unknown"
print("consolidate SF column...")
counts\$Super_family<-paste(counts\$Class, counts\$Super_family, sep = "/")
print("renames NA in SF columns...")
counts\$Super_family[counts\$Super_family == "Unknown/NA"]<-"Unknown"
print(counts[1:20,])
print("compute ecp...")
ecps<-as.data.frame(cbind(as.numeric(counts[,6])/counts\$REPlen, as.numeric(counts[,7])/counts\$REPlen))
print("rename ecp cols...")
print(ecps[1:10,])
print(c(paste(data1, "_ecp", sep = ""), paste(data2, "_ecp", sep = "")))
#names(ecps)<-c(paste(data1, "_ecp", sep = ""), paste(data2, "_ecp", sep = ""))
names(ecps)<-c("ecp_1", "ecp_2")
print("binds to main table...")
counts<-cbind(counts, ecps)


print("exporting table...")
write.table(counts, file="$OUTF/comparison_table.txt", quote = F, row.names = F)


# Plot!

if(ecp_status == FALSE){
print("filtering percent counts...")
counts_t<-as.data.frame(counts[counts[,1] >= pc_T & counts[,2] >= pc_T,])
print("plotting percent...")
plot<-ggplot(na.omit(counts_t), aes(as.numeric(!!ensym(data1)), as.numeric(!!ensym(data2)), col = Class))+
        geom_point()+
        geom_abline(slope = 1, intercept = 0, col = "grey")+
        geom_abline(slope = 10, intercept = 0, col = "red", lty = 3)+
        geom_abline(slope = 0.1, intercept = 0, col = "red", lty = 3)+
        geom_abline(slope = 2, intercept = 0, col = "gold", lty = 2)+
        geom_abline(slope = 0.5, intercept = , col = "gold", lty = 2)+
        scale_y_continuous(trans='log10', breaks = c(0, 0.001, 0.01, 0.1, 1, 10))+
        scale_x_continuous(trans='log10', breaks = c(0, 0.001, 0.01, 0.1, 1, 10))+
        xlab(paste(data1, " (%)", sep = ""))+
        ylab(paste(data2, " (%)", sep = ""))+
        annotation_logticks()+
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
print("export pc plot...")
ggsave(
  paste(data1, data2, "shared_families_percent.pdf", sep = "_"),
  plot = plot,
  device = "pdf",
  path = out,
  scale = 1,
  width = 2600,
  height = 2000,
  units = "px"
)
} else {
print("filtering ecp counts...")
#counts_t<-as.data.frame(counts[counts[,10] >= ecp_T & counts[,11] >= ecp_T,])
counts_t<-counts[counts\$ecp_1 >= ecp_T & counts\$ecp_2 >= ecp_T,]
print("plotting ecp...")
#ggplot(na.omit(counts_t), aes(as.numeric(paste(data1, "_ecp", sep = "")), as.numeric(paste(data2, "_ecp", sep = "")), col = Class))+
plot<-ggplot(na.omit(counts_t), aes(ecp_1, ecp_2, col = Class))+
        geom_point()+
        geom_abline(slope = 1, intercept = 0, col = "grey")+
        geom_abline(slope = 10, intercept = 0, col = "red", lty = 3)+
        geom_abline(slope = 0.1, intercept = 0, col = "red", lty = 3)+
        geom_abline(slope = 2, intercept = 0, col = "gold", lty = 2)+
        geom_abline(slope = 0.5, intercept = , col = "gold", lty = 2)+
        #scale_y_continuous(trans='log10')+ #, breaks = c(0, 0.001, 0.01, 0.1, 1, 10))+
        #scale_x_continuous(trans='log10')+ #, breaks = c(0, 0.001, 0.01, 0.1, 1, 10))+
        xlab(paste(data1, " (equivalent copy)", sep = ""))+
        ylab(paste(data2, " (equivalent copy)", sep = ""))+
        #annotation_logticks()+
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank())
print("export ecp plot...")
ggsave(
 paste(data1, data2, "shared_families_ecp.pdf", sep = "_"),
 plot = plot,
 device = "pdf",
 path = out,
 scale = 1,
 width = 2600,
 height = 2000,
 units = "px"
)
}

SCRIPT
