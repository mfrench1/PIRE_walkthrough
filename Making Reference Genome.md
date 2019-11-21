I'm assuming that you are in the top level of your project directory. You should see your seqeuences ending with `F.fq.gz` and `R.fq.gz` as well as the following directories: `fastq`, `mkREF`, `mkBAM`, `dDocentHPC`, etc...    In the `mkREF` dir are trimmed `fq.gz` and `fq.gz` files that have been stringently filtered by length.  

1. Update the dDocentHPC repository by navigating to it, and pulling down updates which I have successfully verified to work correctly in RPE mode to assemble reference genomes.
  ```bash
  cd dDocentHPC
  git pull
  cd ..
  ```
  *Note* there are 5 primary differences between dDocent's assembly and dDocentHPC mkREF with respect to our single digest data.  

  * First, dDocentHPC removes all singleton sequences within an individual which is computationally advantageous and these sequences will be removed from consideration by the "cutoffs" anyway.  

  * Second, dDocentHPC only assembles contigs that fall between the 5th and 95th percentiles with respect to the number of sequences contained in a precluster. dDocent is set to assemble preclusters with 2 to 10,000 seqs: `rainbow merge -r 2 -R10000`.  For comparison, on _Siganus spinus_ with 15.15 cutoffs, dDocentHPC had the following: `rainbow merge -r 1047 -R 4682` .  I have noticed a relationship between `-R` and the amount of ram being used, which prompted me to implement this change which also prevents us from considering contigs with improbably small or large numbers of sequences. Given that we have to select a smaller set of viable sequences for probes, it seemed ok to me, but I'm open to counter points.

  * Third, I disabled the adapter filter in the middle of dDocentHPC mkREF because (1) it was erroneously generating a few very short sequences that were causing fatal errors downstream and (2) dDocentHPC already filtered for adapters, but dDocent may not have at that point.  

  * Fourth, for speedy reruns at different cutoff values, I have coded dDocentHPC mkREF so that it does not overwrite existing files.  If a file exists and has non-zero size, it will be used.  This does open the possibility for incomplete files to cause failures.  I have tried to improve the feedback the dDocentHPC provides to be able to identify the offending file, and eventually I will add an "overwrite" option.  But for now, if you can't figure out why a run fails, delete all of the files created by dDocentHPC mkREF and start from scratch.

  * Fifth, and this is really minor, but dDocentHPC changed the way that the `rbdiv*out` file is split for running `rainbow merge` in parallel.  Explicitly, dDocentHPC creates a directory called `RBDIV.*` where each `rbdiv*out` file is composed of 1 precluster.  I did this because of exceedingly slow performance and massive ram usage when either running `rainbow merge` in serial, as dDocentHPC used to, or in parallel as `dDocent` does currently with `rebdiv*out` being split into a few hundred files. Obviously different data sets will have varying results here, but this new way results in faster run times, less ram usage, and I've yet to see a node crash like I was with the other strategies.  I recommend, however, against trying to look at the `RBDIV.*` dir directly, which can have over 100,000 files.

2. Copy the `dDocentHPC` files to the `mkREF` directory and move there.
  ```bash
  cp ../dDocentHPC.bash mkREF    #get newest dDocentHPC.bash
  cp config.4.all mkREF                   #copy your config file to build upon the trim settings
  cp dDocentHPC.sbatch mkREF    #copy your sbatch file which has already been customized to your system
  cd mkREF
  ```

3. Update the `config.4.all` file.  Set "Type of reads" to RPE and choose high cutoffs because they will run more quickly and allow you to identify problems more quickly.  Once you get a valid ref genome with the high cutoffs, try some low cutoffs.  For _Siganus spinus_, cutoffs of 20.20 resulted in nearly 1000 contigs and finished within several hours.  15.15 finished within 24 hours.  10.10 and 5.5 are still running. *Note the two lines that I added to the config file, you can copy and paste them into your config file.  They are also in the updated dDocentHPC repo*

  ```bash
----------mkREF: Settings for de novo assembly of the reference genome--------------------------------------------
PE		Type of reads for assembly (PE, SE, OL, RPE)					PE=ddRAD & ezRAD pairedend, non-overlapping reads; SE=singleend reads; OL=ddRAD & ezRAD overlapping reads, miseq; RPE=oregonRAD, restriction site + random shear
0.9		cdhit Clustering_Similarity_Pct (0-1)							Use cdhit to cluster and collapse uniq reads by similarity threshold
10		Cutoff1 (integer)												Use unique reads that have at least this much coverage for making the reference	genome
10		Cutoff2 (integer)												Use unique reads that occur in at least this many individuals for making the reference genome
0.05	rainbow merge -r <percentile> (decimal) 						Percentile-based minimum number of seqs to assemble in a precluster
0.95	rainbow merge -R <percentile> (decimal)							Percentile-based maximum number of seqs to assemble in a precluster
------------------------------------------------------------------------------------------------------------------
  ```

4. Update the `dDocentHPC.sbatch` file to have the following line to run `mkREF`

  ```bash
  dDocentHPC.bash mkREF config.4.all
  ```

5. Assemble reference genome

  ```bash
  sbatch dDocentHPC mkREF
  ```

6. Troubleshooting.  Evaluate the `*out` file generated by your super computer which contains all the feedback from dDocentHPC.  If you are just running in bash, you could redirect the feedback into an out file.

7. Post processing.  I've noticed that about 10% of my contigs have dinucleotide microsats.  To remove them:

```bash
echo -e "CACACACACACA\nTATATATATATA\nGAGAGAGAGAGA\nCTCTCTCTCTCT\nGTGTGTGTGTGT\nGCGCGCGCGCGC" > microsat.motifs
grep -vB1 -f  microsat.motifs reference.15.15.fasta > reference.noSTR.15.15.fasta
```
  * If we think this is a good idea, perhaps somebody can track down a more complete list of motifs to remove and I can figure out how to make `agrep` do this. 
  
  I'm running into an issue with running out of ram (I have 256 gb) when allowing rainbow merge to run in parallel (the default of dDocent). The nature of this single digest data defeats many of the tricks that Puritz used ( and I took even farther) to make this reference assembly work so efficiently. I am using 95GB of ram with 1 process right now, so I'm trying some different things to make this work and finish more quickly.

I am reducing the the -N setting based upon a histogram I made from the rbdiv*out file because I was crashing the node with 40 threads and -N 10000 -R 10000. Here's how to make the histogram for your species (note that you may have different numbers or no numbers in the file name below):

cut -f5 rbdiv.16.16.out | uniq -c | tr -s " " "\t" | sed 's/^\t//g' > rbdiv.16.16.seqsPERprecluster.tsv
  
