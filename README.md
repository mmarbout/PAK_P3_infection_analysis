# PAK_P3_infection_analysis

this set of scripts is dedicated to the analysis of PAK_P3 infection cycle and refered to the publication:
"3D genomic reveal specific patterns of phage genome dynamics and host interactions during lytic infection of Pseudomonas aeruginosa."

## dataset

all the HiC reads can be downloaded on SRA using the following Bioproject reference: PRJNA XXX
The RNA seq reads from Chevalereau et al. 2016 can be downloaded on the NCBI GEO portal (GSE76513).

## dependencies

- Hicstuff (https://github.com/koszullab/hicstuff)
- cooler (https://github.com/open2c/cooler)
- R packages (https://github.com/js2264/HiCExperiment)
- bowtie2

## home made scripts

different scripts need to be downloaded from the github repository.

## analysis

### contact map construction (mcool files)

different argument need to be set (you will have to do that for the different libraries)

genome=PAK_PAKP3
project=MM37
reads_for=MM37_nvq_R1.fq.gz
reads_rev=MM37_nvq_R2.fq.gz
out_reads=MM37_digest
cpu=64

first, pre-digest the reads using the cutsite module from hicstuff.

```sh
hicstuff cutsite -t "$cpu -1 "$reads_for" -2 "$reads_rev" -e "$enzyme" -p "$out_reads"
```
this command will create two fastq files named: MM37_digest_R1.fq.gz and MM37_digest_R2.fq.gz

Then, prepare the bowtie index for the alignment of the reads.

```sh
bowtie2-build "$genome" "$project"/index/idx
```
Launch the pipeline hicstuff to generate the mcool files 

```sh
hicstuff pipeline -t "$cpu" -D -d -f --plot -F --binning 500 --enzyme DpnII,HinfI -o "$project" -g "$project"/index/idx	"$out_reads"_R1.fq.gz "$out_reads"_R2.fq.gz
```
at this step, you should obtain a mcool file: MM37.mcool

You have to obtain all the mcool files for the different libraries.
Several libraries have been sequenced on different lanes of sequencing apparatus (i.e. MM142, MM143 ...) and will need to be merged using cooler.

```sh
cooler merge lib_merge_500bp.cool lib1.mcool::/resolutions/500 lib2.mcool::/resolutions/500 lib3.mcool::/resolutions/500
cooler zoomify -r 1000N --balance -o lib_merge.mcool lib_merge_500bp.cool
```

### RNAseq track generation (bigwig files)


### contact map generation for host and phage

contact map generation for the P. aeruginosa PAK (you can generate contact map at different resolutions)

```sh
plot_matrices_zoom.R lib.mcool lib_PAK_5kb.pdf PAK:1-6000000 5000 
```

contact map generation for the PAK_P3 phage 

```sh
plot_matrices_zoom.R lib.mcool lib_PAK_5kb.pdf PAK_P3:1-88000 500 
```

### Directionnal index


### pileup HiC - RNA

### RNA track plots

### HiC and RNA correlation


### 3D model

### contact map comparison

### 4C-plot

### 4C and RNA correlation







