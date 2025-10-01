# PAK_P3_infection_analysis

this set of scripts is dedicated to the analysis of PAK_P3 infection cycle and refered to the publication:
"3D genomic reveals phage genome dynamics and host interactions during lytic infection of Pseudomonas aeruginosa."

## dataset

all the HiC reads can be downloaded on SRA using the following Bioproject reference: PRJNA1331554
The RNA seq reads from Chevalereau et al. 2016 can be downloaded on the NCBI GEO portal (GSE76513).
The FastA reference file can be downloaded on the github.

if you want to skip the generation of raw data files, you can find them at the following link : 
- faire un repo Zenodo pour les fichiers mcool & co

## dependencies

- Hicstuff (https://github.com/koszullab/hicstuff)
- cooler (https://github.com/open2c/cooler)
- R packages (https://github.com/js2264/HiCExperiment)
- bowtie2
- tidycoverage (https://github.com/js2264/tidyCoverage)

## home made scripts

different scripts need to be downloaded from the github repository and put in a dedicated directory named [scripts]

## Data generation

the first steps will be to generate the HiC mcool files and the RNA bigwig files.

### contact map construction (mcool files)

different argument need to be set (you will have to do that for the different libraries)

- genome=PAK_PAKP3.fa
- project=MM37
- reads_for=MM37_nvq_R1.fq.gz
- reads_rev=MM37_nvq_R2.fq.gz
- out_reads=MM37_digest
- cpu=64

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
at this step, you should obtain a mcool file: MM37.mcool as well as other mandatory files to pursue the analysis like the distance law file: distance_law.txt

You have to obtain all the output files for the different libraries.
Several libraries have been sequenced on different lanes of sequencing apparatus (i.e. MM142, MM143 ...) and will need to be merged using cooler.

```sh
cooler merge lib_merge_500bp.cool lib1.mcool::/resolutions/500 lib2.mcool::/resolutions/500 lib3.mcool::/resolutions/500
cooler zoomify -r 1000N --balance -o lib_merge.mcool lib_merge_500bp.cool
```

### RNAseq track generation (bigwig files)

```sh
tinyMapper.sh -m RNA -s "$work_dir"/fastq/RNA/"$project" -o "$work_dir"/RNA_track/"$project" -g "$work_dir"/fasta/deinoc_V2	--threads 24
```

## analysis

### contact map generation for host and phage

contact map generation for the P. aeruginosa PAK (you can generate contact map at different resolutions)

```sh
plot_matrices_zoom.R lib.mcool lib_PAK_5kb.pdf LR657304.1:1-6395872 5000 
```

contact map generation for the PAK_P3 phage 

```sh
plot_matrices_zoom.R lib.mcool lib_PAKP3_500b.pdf NC_022970.1:1-88097 500 
```

you can obviously do a for loop to generate all the contact map for the different libraries and at different resolutions

```sh
for lib in MM53 MM54 MM55 MM56 MM57 MM58
do
  for resolution in 1000 2000 5000
  do
    plot_matrices_zoom.R "$lib".mcool "$lib"_PAK_5kb.pdf LR657304.1:1-88097 "$resolution"
  done
  for resolution in 500 1000
  do
     plot_matrices_zoom.R "$lib".mcool "$lib"_PAKP3_5kb.pdf NC_022970.13:1-88097 "$resolution"
  done
done
```


### Directionnal index

DI analysis allow to detect CIDs in bacterial chromosomes (i.e. Chromosome Interacting Domains)

```sh
directionnal_index.R MM37.mcool LR657304.1:1-6395872 10000 10 Fig1+supp/CID_MM37.pdf 
```
### HiC and RNA correlation

it is also possible to diretly look at the correlation between HiC signal and transcription signal

```sh
HiC_RNA_correlation.R MM37_filter.mcool 1000 LR657304.1:1-6395872 T0_rep1_unstranded.bw 1000 MM37_HiC_RNA_correlation.pdf
```



### pileup HiC - RNA

### RNA track plots


### distance law plots




### 3D model

### contact map comparison

### 4C-plot

### 4C and RNA correlation







