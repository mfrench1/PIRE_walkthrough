[dDocentHPC](https://github.com/cbirdlab/dDocentHPC) will be used to trim the files.  

dDocentHPC is like dDocent but it has been modified to behave more like a unix program (disabled interactivity). The "config" file settings have also be reworked. All settings are controlled by the config file and the script dDocentHPC.bash must be provided the config file as an argument at the command line when running.  There is also a dDocent.sbatch file to run the script on a SLURM cluster. For trimming, this version has been modified to work like earlier versions of dDocent that trimmed differently for reference creation and mapping. This allows for more stringent filtering to be applied for reference creation and more leniently filtering for mapping.

In terms of software, dDocentHPC and dDocent have the same dependencies (running on the software required by dDocent 2.7.8)

Assuming that you are starting in your working directory and the fq.gz files have been prepared and named as described below in the "Demultiplexing Files" section.

Clone GitHub repository 
```
git clone https://github.com/cbirdlab/dDocentHPC.git
```

Copy files from repo to _Species Dir_
```
cp dDocentHPC/dDocentHPC* .
cp config.4.all .
```

Edit the `config.4.all` file.  You should customize the first two values (#processors, MaxMemory) to your system's specifications.  The third setting controls read length to be used for reference creation.  In the example below, 140 was used because the `F.fq.gz` files are 141 bp long.  If there is more than 1 base clipped, then the whole paired-end read will be tossed.  There can be up to 11 bp clipped from the `R.fq.gz` reads which are 151 bp long.  The other settings have worked well for us and probably don't need to be changed.  If we lose a substantial proportion of reads with the stringent 140 bp threshold, then we can go back and relax it until a balance is struck.  The threshold is more lenient for mapping with reads below 75 bp being discarded.  
```
Settings for dDocentHPC
These default settings assume ddRAD, no overlapping 151 bp reads

40              Number of Processors (Auto, 1, 2, 3, ..., n threads) cbirdq=40 normal=20
230G            Maximum Memory (1G,2G,..., 256G)  G=gigabytes

----------trimFQ: Settings for Trimming FASTQ Files---------------------------------------------------------------
140             trimmomatic MINLEN (integer, mkREF only)                        Drop the read if it is below a specified length. Set to the length of the Read1 reads.
75              trimmomatic MINLEN (integer, mkBAM only)                        Drop the read if it is below a specified length. Set to the minimum frag length you want$
20              trimmomatic LEADING:<quality> (integer, mkBAM only)             Specifies the minimum quality required to keep a base.
15              trimmomatic TRAILING:<quality> (integer, mkREF only)            Specifies the minimum quality required to keep a base.
20              trimmomatic TRAILING:<quality> (integer, mkBAM only)            Specifies the minimum quality required to keep a base.
2               trimmomatic ILLUMINACLIP:<seed mismatches> (integer)            specifies the maximum mismatch count which will still allow a full match to be performed
30              trimmomatic ILLUMINACLIP:<palindrome clip thresh> (integer)     specifies how accurate the match between the two 'adapter ligated' reads must be for PE $
10              trimmomatic ILLUMINACLIP:<simple clip thresh> (integer)         specifies how accurate the match between any adapter etc. sequence must be against a rea$
20              trimmomatic SLIDINGWINDOW:<windowSize> (integer)                specifies the number of bases to average across
20              trimmomatic SLIDINGWINDOW:<windowQuality> (integer)             specifies the average quality required.
0               trimmomatic HEADCROP:<length> (integer, only Read1 for ezRAD)   The number of bases to remove from the start of the read. 0 for ddRAD, 5 for ezRAD
no              FixStacks (yes,no)                                                                                      Demultiplexing with stacks introduces anomolies.$------------------------------------------------------------------------------------------------------------------
```

Edit the `dDocentHPC.sbatch` file. One job is to trim the reads for reference creation and one is to trim the reads for mapping. It goes faster if you run them at the same time rather than consecutively so just modify  `dDocentHPC.sbatch` for each new run.  Make sure to comment out most of the lines, except those listed here:
```
#!/bin/bash

#SBATCH --job-name=trimFQ
#SBATCH -o trimFQref-%j.out
#SBATCH --time=48:00:00
#SBATCH -p cbirdq
#SBATCH --nodes=1

#enable lmod
module load ddocent/2.7.8

bash dDocentHPC.bash trimFQref config.4.all
#bash dDocentHPC.bash trimFQmap config.4.all
```

On the ODU HPC, you need a ` -l` after the shebang!, ` #SBATCH -p main`, delete `#SBATCH --nodes=1`, add `#SBATCH --ntasks=32`, and uncomment `enable lmod` before module load, if I remember correctly.

Run the script.  It will create two dirs: `mkREF` and `mkBAM`.  The files trimmed for reference creation will go into `mkREF` and those for mapping will go into `mkBAM`.  A `BAM` file, binary alignment map, is the file that results from mapping, thus the dir name.  
```
sbatch dDocentHPC.sbatch
```

Review the  *.out files for errors.

Make sure that all of the files were trimmed.  I you get a "java ran out of memory" error, dial back the number of threads. 

Evaluate the proportion of reads lost in the trimming.  If a lot were lost, go back, lower the `140` setting and try again.
