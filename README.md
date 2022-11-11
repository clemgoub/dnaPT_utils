# dnaPT_utils

dnaPT_utils is a collection of scripts to perform multiple downstream analyses with [dnaPipeTE 1.3](https://github.com/clemgoub/dnaPipeTE)
I recommend to use this tools to replace the graphs automatically produced by the pipeline.

- [`dnaPT_charts.sh`](#dnapt_chartssh): plots the relative repeat proportions (2 piecharts) and the quantification of each assembled repeat (barplot).
- [`dnaPT_landscapes.sh`](#dnapt_landscapessh): plots an histogram representing the blastn divergence mesured between each read and the assembled repeats.
- [`dnaPT_compare.sh`](#dnapt_comparesh): perform a comparative analysis of the relative abundance of shared repeat families between two datasets.

Contributors:
- T. Mason Linscott

## dnaPT_charts.sh

This script processes an output folder of dnaPipeTE v1.3 and produces 3 graphs:
- Piechart 1 (left): relative proportions of the different repeat categories.
- Piechart 2 (right): proportion of the different repeat categories relative to the total genome.
- Barplot (bottom): genomic proportion of each dnaPipeTE contig and the associated classification.
	- A threshold `-t` (% genome) can be set such as contigs representing < t% of the genome will be groupped together.

![charts](./figures/charts.png)
> *You will note that the piecharts may provide slightly different values than originaly. Indeed, the piechart directly available after a dnaPipeTE run is based on the counts of 2 successive blastn again (1) the classified repeats in `Trinity.fasta` and (2) the unclassified ones (these counts can be found in `Counts.txt`). This was originaly designed to increase sensitivity at the subclass level. However, it turns out to be confusing without providing  substantial improvement. With this script, the piecharts are now directly made from the results of one single blastn of the reads against the whole `Trinity.fasta` file, which is also used to create the detailed count table `reads_per_components_and_annotation` as well as the landscape analysis.*

### Dependencies:
- `R` + packages `ggplot2`, `tidyverse` and `cowplot` (https://www.r-bloggers.com/2010/11/installing-r-packages/)

### Usage: 

```
./dnaPT_charts.sh -I <dataset_directory> [options]


mendatory arguments:
 -I, --input-dir              dnaPipeTE output directory

options:
 -p, --prefix                 prefix to append to the output filename: "<prefix>_charts.pdf"
 -t, --percent_threshold      barplot: repeats < -t % are groupped together as a single category (default 0.001%)
 -o, --output                 output folder (path); default: dnaPipeTE output directory
 -y, --ylim                   Max value for the y axis (genome %) [0-100]
 -h, --help                   Prints this message and exit
 ```

## dnaPT_landscapes.sh

This script perform at "TE landscape" analysis, i.e., it plots an histogram of the blastn divergence between raw reads (TE copies in the genomes) and their consensus sequences assembled in the file `Trinity.fasta`. The script plots only putative TE sequences among the subclasses "LINE", "SINE", "LTR", "DNA", "RC" and "Unknown" (a.k.a. "NA"). 

![landscapes](./figures/landscapes.png)

### Dependencies:
- R + packages `ggplot2` and `tidyr`

### Usage: 

```
./dnaPT_landscape.sh -I <dataset_directory> [options]

mendatory arguments:
 -I, --input-dir              dnaPipeTE output directory

options:
 -p, --prefix                 prefix to append to the output filename: "<prefix>_landscapes.pdf"
 -o, --output                 output folder (path); default: dnaPipeTE output directory
 -S, --superfamily            Plot with superfamily information (instead of subclass)
 -y, --ylim                   Max value for the y axis (genome %) [0-100]
 -U, --no-unknown             Remove unclassified repeats
 -h, --help                   Prints this message and exit
```

## dnaPT_compare.sh

This script measures the relative abundance of shared TE families between two datasets analyzed with dnaPipeTE. 
   
1 - Shared families are identified by clustering together the repeat sequences ("`Trinity.fasta`") of each dataset. The clustering is performed with `cd-hit-est` with the parameters described in Goubert et al, 2021 (Mobile DNA, in press) and approximate the "80-80-80" rule. 
   
2 - Shared families are identified by selecting clusters where sequences from both dataset A and B are present. For each cluster, the counts (in bp) and genome % are summed per dataset to obain a quantification of each shared family. The classification of a shared repeat is taken from the representative sequence of each cluster, and correspond to the longest sequence in the cluster. It can either come from dataset A or B. 
   
3 - The abundances of each shared family, either in % genome or equivalent copy, are then plotted with R/ggplot2.

![compare](./figures/compare.png)
>*I recommend to use caution while interpreting results with low quantities, typically < 0.01% or < 1 equivalent copy. Thresholds options are available (-p / -e) to filter the plotted data.*

### Dependencies:

- CD-HIT (https://github.com/weizhongli/cdhit)
- R + packages `ggplot2`, `scales`, `tidyr` and `reshape2` 

### Usage: 

```
./dnaPT_compare.sh -A <dataset_A_directory> -a <prefix_A> -B <dataset_B_directory> -b <prefix_B> -o <output_folder> [options]

mendatory arguments:
 -A, --dir_A                  dnaPipeTE output directory for dataset A (path)
 -B, --dir_B                  dnaPipeTE output directory for dataset B (path)
 -a, --pref_a                 prefix for dataset A (string)
 -b, --pref_b                 prefix for dataset B (string)
 -o, --output                 output folder (path)

options:
 -T, --te_only                Only plot repeats of the sub-classes "LINE", "SINE", "LTR", "DNA", "RC" and "Unknown"
 -S, --superfamily            Plot with superfamily information (instead of Class)
 -p, --percent_threshold      min. percent genome to plot (default = 0) / not used if -E/-ecp selected
 -E, --ecp                    perform comparison in equivalent copy (dataset counts in bp / representative sequence size)                              
 -e, --ecp_threshold          min. equivalent copy/ies to plot (default = 0)

 -h, --help                   Prints this message and exit
```
> *Caution: the length used for normalization in "equivalent copy" mode (`-E/--ecp`) is based on the length of the representative sequence of each cluster, which originate from either one of the two dataset compared. It is thus assumed that the consensus length is the same in each species/sample, which is not necessarily true.*

### Output files
In addition to the graphs exported in pdf, this script also produce the following files that may be useful for the user (other files are intermediates). A and B will be replaced by the dataset prefix:
- `A_B`: fasta file with only the representative sequences of each cluster (produced by `cd-hit-est`)
- `A_B.clstr` and `A_B.bak.clstr`: clustering output of `cd-hit-est` (first is the main output, the second is a table version of it)
- `A_B_comparison_table.txt`: processed data table used for plotting in R, unfiltered (i.e. is produced before applying the filtering options `-T`, `-p` or `-e`). The two first columns, named A and B according to the datasets prefixes, report the % genome of each shared repeats. The two last columns, ecp_1 and ecp_2, report the "equivalent copy" for A and B respectively.
- `A_B_R.tsv`: raw data table (used as input for the R code)
- `A_B_dnaPipeTE_contigs.fasta`: concatenated `Trinity.fasta` files from each dataset, with dataset prefix added. Input for `cd-hit-est`.

### Known issues:

- A common source of error is if one of more annotations are missing in the file `colors.land`:

This error causes the scripts `dnaPT_charts.sh` and `dnaPT_landscape.sh <with -S option>` to fail with the error: 

```
Error in col.bars[i] <- cols$V2[grep(pattern = paste("^", levels(as.factor(rca_t$t_subc))[i],  :
  replacement has length zero
Execution halted
```
> See also issue #4

The solution is to add the missing annotation to `colors.land`. The file is present alongside the `dnaPT_utils` scripts and look like so (extract):
```
DNA/TcMar-Tigger    "#FF936C"
DNA/TcMar-Tigger?   "#FF936C"
DNA/Zator   "#FF9F79"
DNA/Zator?  "#FF9F79"
DNA/Zisupton    "#FFCCCC"
DNA/Zisupton?   "#FFCCCC"
LINE    "royalblue"
LINE?   "#251792"
LINE/Ambal  "#37B9F0"
LINE/CR1    "#483AA2"
LINE/CR1?   "#483AA2"
```

To add your annotation, follow this pattern:

- 1 If the annotation is of the for `class/superfamily` (or similar):
    - 1.1 add one line `class/superfamily<tab>color` **AND**
    - 1.2 add on line with only the class: `class<tab>color`.
- 2 If the annotation is only `class` (no `/`):
    - 2.1 only add the line `class<tab>color`

`color` needs to be in hexadecimal form with double quotes (e.g. `"#AD49F2"`) see: https://htmlcolorcodes.com/ **OR** a R color (e.g. `"royalblue"`) see: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.

To find the missing annotations, the script will output the latest annotations before the error like so:

```
[1] "Recognized classes:"
[1] "^DNA$"
[1] "^LINE$"
[1] "^Low_complexity$"
[1] "^LTR$"
[1] "^PLE$"
Error in col.bars[i] <- cols$V2[grep(pattern = paste("^", levels(as.factor(rca_t$t_subc))[i],  :
  replacement has length zero
Execution halted
```
> In this case, the first annotation that was missing was `PLE` (it is now supported)
