I'm assuming that your are in the top level of your project directory. You should see your seqeuences ending with F.fq.gz and R.fq.gz as well as the following directories: fastq, mkREF, mkBAM, dDocentHPC, etc... In the mkBAM dir are trimmed fq.gz and fq.gz files that have been filtered by length (min of 75 bp).

*Eric's note: Chris has recommended running the "Barcode-Restriction Site Motif" filter before mapping. Also, in later steps we have been changing the names of some files to reflect the step in which files are.  If you do run the "Barcode-Restriction Site Motif" before mapping, make sure not to change the name of the reference file (reference*.fasta) or the script won't run correctly.

1. Update the dDocentHPC repository by navigating to it, and pulling down updates which I have successfully verified to work correctly in RPE mode to assemble reference genomes.
```bash
cd dDocentHPC
git pull
cd ..
```

2. Copy necessary files to the mkBAM directory.
```bash
cp mkREF/config.4.all mkBAM   #make sure you have the latest version with rainbow merge settings file 09-29-2019
cp dDocentHPC/dDocentHPC.bash mkBAM
cp mkREF/dDocentHPC.sbatch mkBAM
cp mkREF/reference*fasta mkBAM
```

I have modified dDocent's mapping heavily in the past, primarily with regards to settings.  I'm presently unsure how those changes will affect RPE (single digest) data.  I know that they create beautiful ddRAD data.  What you should know is that dDocentHPC makes the BAM files then separately filters them `dDocentHPC fltrBAM`,  which is different than dDocent which makes the BAM files and filters them together and you easily can't adjust  the filtering.  Further, the default dDocentHPC mkBAM and fltrBAM settings are different than dDocent.

3. Move the mkBAM dir and update the config.4.all and sbatch files accordingly.  While the cutoff settings occur in the assembly portion of the settings, they are used throughout, so make sure they match your files.
```bash
#code to run mapping in sbatch file
bash dDocentHPC.bash mkBAM config.4.all
```

4. Run the mapping
```bash
sbatch dDocentHPC.sbatch
```
