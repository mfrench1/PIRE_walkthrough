I'm assuming that you're starting in your mkBAM dir

1. Create a mkVCF directory and move the non-"Pool" `*RG.bam` files (fq.gz and bam) to it
```bash
mkdir ../mkVCF
ls *RG.bam | grep -v '_Pool-' | parallel --no-notice "mv {} ../mkVCF"
cd ../mkVCF
```
1.2. Copy reference genome to mkVCF

2. Pull down the latest version of dDocentHPC.  I made several changes on 09-29-2019, and it will be important to get both the newest config.4.all and the dDocentHPC.bash files.

3. Edit config file.  I made it so that all positions will be typed and it will mostly only look at read1 when genotyping.  That's not always the case right now, but it is mostly the case.  I'll see if I can improve the script as necessary or we can cull positions from the *.vcf.  I have made it so that SNPs, insertions, and deletions are called while MNPs (multinucleotide polymorphisms) are called as SNPs.

```bash
----------mkVCF: Settings for variant calling/ genotyping---------------------------------------------------------
no              freebayes -J --pooled-discrete (yes or no)                                              If yes, a pool of individuals is assumed to be the statistical unit of observation.
no              freebayes -A --cnv-map (filename.bed or no)                                             If the pools have different numbers of individuals, then you should provide a copy number variation (cnv) *.bed file with $
2               freebayes -p --ploidy (integer)                                                                 Whether pooled or not, if no cnv-map file is provided, then what is the ploidy of the samples? for pools, this num$
no              freebayes -r --region (filename.bed or no)                                              Limit analysis to specified region.  Bed file format: <chrom>:<start_position>-<end_position>
150             only genotype read 1 (integer)                                                                  Limit analysis to only Read 1 positions, integer is maximum Read1 bp position
0               freebayes -n --use-best-n-alleles (integer)                                             reduce the number of alleles considered to n, zero means all, set to 2 or more if you run out of memory
30              freebayes -m --min-mapping-quality (integer)
20              freebayes -q --min-base-quality (integer)
-1              freebayes -E --haplotype-length (-1, 3, or integer)                     Set to -1 to avoid multi nucleotide polymorphisms and force calling MNPs as SNPs.  Can be set up to half the read length, or more.
0               freebayes    --min-repeat-entropy (0, 1, or integer)                    Set to 0 to avoid multi nucleotide polymorphisms and force calling MNPs as SNPs. To detect interrupted repeats, build across sequence unti$
10              freebayes    --min-coverage (integer)                                           Require at least this coverage to process a site
0.25    freebayes -F --min-alternate-fraction (decimal 0-1)                     There must be at least 1 individual with this fraction of alt reads to evaluate the position. If your individuals are barcoded, then use 0.2. If y$
3               freebayes -C --min-alternate-count (integer)                                    Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the pos$
3               freebayes -G --min-alternate-total (integer)                                    Require at least this count of observations supporting an alternate allele within the total population in order to use the allele $
0.33    freebayes -z --read-max-mismatch-fraction (decimal 0-1)                 Exclude reads with more than N [0,1] fraction of mismatches where each mismatch has base quality >= mismatch-base-quality-threshold default: 1.0
20              freebayes -Q --mismatch-base-quality-threshold (integer)                Count mismatches toward --read-mismatch-limit if the base quality of the mismatch is >= Q.  default: 10
50              freebayes -U --read-mismatch-limit (integer)                    Exclude reads with more than N mismatches where each mismatch has base quality >= mismatch-base-quality-threshold. default: ~unbounded
20              freebayes ~3 ~~min-alternate-qsum (integer)                                             This value is the mean base quality score for alternate reads and will be multiplied by -C to set -3. Description of -3: R$
50              freebayes -$ --read-snp-limit (integer)                         Exclude reads with more than N base mismatches, ignoring gaps with quality >= mismatch-base-quality-threshold. default: ~unbounded
20              freebayes -e --read-indel-limit (integer)                                               Exclude reads with more than N separate gaps. default: ~unbounded
no              freebayes -w --hwe-priors-off (no|yes)                                                  Disable estimation of the probability of the combination arising under HWE given the allele frequency as estimated by obse$
no              freebayes -V --binomial-obs-priors-off (no|yes)                                 Disable incorporation of prior expectations about observations. Uses read placement probability, strand balance probability, and r$
no              freebayes -a --allele-balance-priors-off (no|yes)                               Disable use of aggregate probability of observation balance between alleles as a component of the priors.
no              freebayes --no-partial-observations (no|yes)                                    Exclude observations which do not fully span the dynamically-determined detection window.  (default, use all observations, dividin$
yes             freebayes    --report-monomorphic (no|yes)                                              Report even loci which appear to be monomorphic, and report allconsidered alleles, even those which are not in called geno$
------------------------------------------------------------------------------------------------------------------

```

4. Update the `dDocentHPC.sbatch` file with the `dDocentHPC.bash mkVCF` command

5. run
