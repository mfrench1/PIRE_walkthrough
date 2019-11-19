Files should be demultiplexed using the following parallel demultiplex script.  https://github.com/cbirdlab/RAD_demultiplex

# *RAD_demultiplex

Scripts to demultiplex internally barcoded fastq files on a SLURM scheduled HPC cluster in parallel.  With some mild hacking, this will run on a workstation in bash.  Each pair of fastq files to be demultiplexed can be run on a different thread if multiple files are specified in the `DMXfiles` and `FQfiles` variables, see *To Run* below. Demultiplexed fastq files are saved to directories in the `pwd`.

ddRAD_demultiplex.sbatch is for double-digest RAD data

newRAD_demultiplex.sbatch is for single-digest RAD data that has 1 ligated barcode followed by ezRAD sequencing (aka newRAD or bestRAD)
* Note that `process_radtags` which is harnessed by this script cannot handle
  * The remaining restriction site before the barcode, for SbfI, it's `GG`, so it should be added to the barcode in the demultiplex decode files
  * *_Indels at the beginning of the sequence reads, which are quite common. I have to code a solution to this using `agrep ` _*

## To Run

* Prepare your files.  There should be 1 demultiplex decode file per pair of fastq files (assuming paired end sequencing, R1 & R2), and each should be formatted with the first column being the barcodes, a tab, and the second column being the base name of the resulting demultiplexed sequences:
  ```
  GGAAGCCGGT      PIRE2019-Ssp-C-Gub_096-Plate1Pool6Seq1-2G-L4
  GGCGATGCTC      PIRE2019-Ssp-C-Gub_068-Plate1Pool6Seq1-2G-L4
  ```
  The base names of the demultiplex decode files should match those of the fq files they refer to, and the code assumes that the name of the   demultiplex decode files ends with `_demultiplex.txt`

* Clone this repo to your computer
  ```
  git clone https://github.com/cbirdlab/RAD_demultiplex.git
  ```

* Copy the appropriate script to the directory where you want the demultiplexed files to be saved.

* No commandline arguments are accepted.  To specify where your data is, edit the following variables which hold the paths to the demultiplex decode and fastq files as well as the READ1 & READ2 file extensions used:
  ```
  #populate variables
  DMXfiles=*_demultiplex.txt
  FQfiles=*_1.fq.gz
  R1Ext=_1.fq.gz
  R2Ext=_2.fq.gz
  ```

* Edit the following `SBATCH` commands to work with your cluster
  ```
  #SBATCH -p normal
  #SBATCH --nodes=1

  #load software tools
  module load stacks
  module load parallel
  ```

* Run the script
  ```
  sbatch ddRAD_demultiplex.sbatch
  ```
  or
  ```
  sbatch newRAD_demultiplex.sbatch
  ```

### Most Common Motifs in 1st 10bp of newRAD
`GGCGATGCTC`   : expected res site - barcode
`GGAAGCCGGT`  : expected res site - barcode
` GCGATGCTC`      : res site - barcode with 1 del
` GAAGCCGGT`     : res site - barcode with 1 del
`  CGATGCTC`         : res site - barcode with 2 del = barcode
`  AAGCCGGT`        : res site - barcode with 2 del = barcode
`GATCGGAAGA `  : i7 or rev-comp i5
`AAGCAGAAGA`   : rev-comp i7
` ATCGGAAGAG`   : i7 or rev-comp i5 with 1 del
`GGACGTAC..`       : biotin side of ligation 1 adapter

1-2 mismatches are common across these motifs

For the first few species, the `BARCODE` is either `CGATGCTC` or `AAGCCGGT`

### Hackish Solution to Recover Good Barcodes and Those with a Leading Deletion

Here are the details on how to hack together most of the reads with legitimate barcodes presuming that, like the _Siganus spinus_ I've looked at, the majority of discarded legitimate barcodes have a single leading deletion.  It would be pretty easy to update this methodology to include barcodes with 2 leading deletions, but it's probably not necessary. It does not, however address that 2 mismatches are fairly common and a new script using `agrep` or `tre-grep` will give us more control there.  Until then, this is what we have.

I'm assuming that you have a _Species Dir_ for this data and have downloaded the `*fq.gz` and `*_demultiplex.txt` files

*Prep _Species Dir_*
* if your barcodes  in `*_demultiplex.txt` don't begin with `GG` , then add `GG`
    e.g., `ls *_demultiplex | parallel --no-notice -j10 -q sed -i 's/^/GG/' {} `
* make a directory called `dmplxGG` in your _Species Dir_ and *copy* your `*_demultiplex.txt` files there
  ```
  mkdir dmplxGG
  cp *_demultiplex.txt
  ```
* change `GG` to `G` in `YourSpeciesDir/*_demultiplex.txt`
    ```
    ls *_demultiplex | parallel --no-notice -j10 -q sed -i 's/^GG/G/' {} 
    ```
* make a directory called `dmplxG` in your _Species Dir_ and *move*  your `*_demultiplex.txt` files there
  ```
  mkdir dmplxG
  mv  *_demultiplex.txt
  ```
* make a directory called `fastq` in your _Species Dir_ and move your `*fq.gz` files there
  ```
  mkdir fastq
  mv *fq.gz fastq
  ```

*_Demultiplexing_*
* Demultiplex in `dmplxGG` dir using [newRAD_demultiplex.sbatch](https://github.com/cbirdlab/RAD_demultiplex)
  ```
  git clone https://github.com/cbirdlab/RAD_demultiplex.git
  cp RAD_demultiplex/newRAD_demultiplex.sbatch dmplxGG
  cd dmplxGG
  nano dmplxGG #modify to work with your system and dirs
  sbatch newRAD_demultiplex.sbatch
  ```
* demultiplex in `dmplxG` dir using [newRAD_demultiplex.sbatch](https://github.com/cbirdlab/RAD_demultiplex)
  ```
  cd ../dmplxG
  cp dmplxGG/newRAD_demultiplex.sbatch .
  sbatch newRAD_demultiplex.sbatch
  ```

*_Combine Demultiplexed Files_*
* Move to _Species Dir_ and Concatenate the `fq.gz` files from `dmplxG` and `dmplxGG`
  ```
  cd ..

  # you have to modify these next 2 lines to match your species file names
  ls dmplxG/demultiplexed_seqs_20190905_PIRE_Ssp-*/*L4.[12].fq.gz > dmplxG.txt
  ls dmplxGG/demultiplexed_seqs_20190905_PIRE_Ssp-*/*L4.[12].fq.gz > dmplxGG.txt
  cat dmplxGG.txt | parallel -k --no-notice "basename " | sed 's/1\.fq/F\.fq/' | sed 's/2\.fq/R\.fq/' > fqnames.txt

  # use parallel and zcat to finish creating the fq.gz files
  parallel --no-notice --link " zcat {1} {2} | gzip > {3} " :::: dmplxGG.txt :::: dmplxG.txt :::: fqnames.txt
  ```
  
  
