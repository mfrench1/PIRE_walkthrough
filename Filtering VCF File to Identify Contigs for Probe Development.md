Our goal is to minimally filter our data, minimize ascertainment bias, and identify contigs to use for probe development.  Here's a brief overview of the filtering steps, in order:

* Fltr30: Eliminate positions from consideration due to low cvg, abnormally high rates of SNPs or indels.  Essentially we snip off the first few bases of each contig, and cut off the bases around position 140 (150bp read - barcode)
  * Reasoning: we don't want to keep the first positions because they almost always have problems, and past ~140 bp, we are getting into lower coverage due to either the unsequenced positions between Read 1 and Read 2 or lower coverage due to the random shearing of the Read 2 side of the DNA.

* Fltr14: Eliminate genotypes from consideration with less than 10x coverage.  We explicitly change every genotype and the associated genotype data in the vcf that had fewer than 10 reads to missing data.
  * Reasoning: since we probably won't really be using genotypes that have less than 10x coverage, this gives us a realistic view of the data for analysis and enables the filtering of individuals that didn't sequence well in the next filter.

* Fltr16: Remove individuals with more than 50% missing data.
  * Reasoning: thus far we have allowed poorly sequenced individuals to move through the pipeline because they don't add much time to the process, but now we need to remove them so we can get good estimates of coverage.

* Fltr31: Remove remaining contigs that are shorter than a given cutoff (around 110-120 bp)
  * Reasoning: I noticed that most contigs represented by a small number of bp have low coverage.  While some of these contigs may be ok, the majority are flawed, so we remove them, leaving us with contigs that have consistently long sequence reads.

* Fltr041: Remove contigs with abnormally low and high mean depth of coverage using percentile cutoffs (0.01 and 0.99) and abnormally high variation in depth of coverage across individuals.  
  * Reasoning: Those with low depth of coverage are not being sequenced reliably and those with high depth of coverage are probably paralogs.  I used the 0.01 and 0.99 because it removes only the most extreme loci and the lower mean cvg cutoff was around 15x, which seems good.  Those with a high degree of variation in cvg among individuals may be due to issues with allelic dropout.

* Fltr32: Remove contigs with abnormally high levels of heterozygosity.
  *Reasoning: If two invariant loci that are different are mapped to the same contig, they will present with every individual being heterozygous for every SNP.  rad_haplotyper does not catch these.  I'm not finding many, so it seems good.
  * Potential improvements: I could make a filter based upon Garrett McKinney's D vs Het methodology for identifying paralogs (but it was easier to run rad_haplotyper)

* Fltr05,07,19: Prepare data for rad_haplotyper by removing positions with too much missing data, then those that are monomorphic.  Then run rad_haplotyper with settings so lenient that it won't filter nearly anything, so that we can use the stats output to manually choose a cutoff for "number of possible paralogs".  Those contigs can then be manually removed from the `*Fltr32*vcf`  with `grep`.
  *Reasoning: we don't want to enrich for paralogous loci with the probes.

* Final step: use final list of contigs to subset the reference genome, yielding the contigs for probe selection.  I strongly suggest that we design probes that target the first 120-250 bp of the contigs.  There's another discussion to be had about how many probes to develop per contig (I'm thinking 2-4).

---

I'm assuming that you're starting in your mkVCF dir

1. Create a `fltrVCF` dir 
```bash
mkdir ../fltrVCF
```

2. Clone the fltrVCF and my fork of the rad_haplotyper repos above the level of your project dir.
```bash
cd ../../
git clone https://github.com/cbirdlab/fltrVCF.git
git clone https://github.com/cbirdlab/rad_haplotyper.git
```

3. Copy a the files necessary for running fltrVCF to the fltrVCF dir inside your project dir
```bash
cp rad_haplotyper/*pl YOURPROJECTDIR/fltrVCF
cp fltrVCF/*pl YOURPROJECTDIR/fltrVCF
cp fltrVCF/fltrVCF* YOURPROJECTDIR/fltrVCF
cp fltrVCF/config* YOURPROJECTDIR/fltrVCF
cp fltrVCF/scripts/*R YOURPROJECTDIR/fltrVCF
cp fltrVCF/scripts/*fltrVCFstats2* YOURPROJECTDIR/fltrVCF
```

4. Set up your config file. Everybody's data is different, so some customization of settings may be necessary.  The goal with the filter 19 settings (rad_haplotyper) is to filter nothing on the first run through
```bash
        fltrVCF -f 30 14 16 31 041 32 05 07 19
        fltrVCF -c 15.15
        fltrVCF -b ../mkVCF/                           #path to *.bam files
        fltrVCF -d ../mkBAM/mapped.15.15.bed
        fltrVCF -v ../mkVCF/TotalRawSNPs.15.15.vcf
        fltrVCF -g ../mkVCF/reference.15.15.fasta
        fltrVCF -p ../mkVCF/popmap.15.15
        fltrVCF -w filter_hwe_by_pop_HPC.pl
        fltrVCF -r rad_haplotyper.pl
        fltrVCF -o PIRE_SiganusSpinus.A
        fltrVCF -t 20                                           #number of threads [1]

        30 custom bash                  8                 #Keep sites after this position (bp)
        30 custom bash                  142             #Keep sites before this position (bp)
        31 custom bash                  113             #Remove contigs with fewer BP
	32 custom bash		0.9		#Keep contigs with lesser porportion of heterozygotes

	041 custom bash				0.05	#Remove contigs with lower mean of mean depth across sites, percentile [0.01]
	041 custom bash				0.95	#Remove contigs with higher mean of mean depth across sites, percentile [0.99]
	041 custom bash				0		#Remove contigs with lower CV of mean depth across sites, percentile [0]
	041 custom bash				0.95	#Remove contigs with higher CV of mean depth across sites, percentile [0.99]
        05 vcftools --max-missing       0.8           #Remove sites with lower proportion of genotypes present [0.5]

	07  vcffilter AC min		0		#Remove sites with equal or lower MINOR allele count [1]

        14 vcftools --minDP             5               #Code genotypes with lesser depth of coverage as NA [5] 

        16 vcftools --missing-indv      0.5     #Remove individuals with more missing data. [0.5]

        19 rad_haplotyper       -d      50              #depth of sampling reads for building haplotypes. [50]
        19 rad_haplotyper       -mp     25              #Remove sites with more paralogous indivduals. Adjust according to sample size. [10]
        19 rad_haplotyper       -u      150              #Remove contigs with more SNPs. Adjust according to sequence length. [30]
        19 rad_haplotyper       -ml     25              #Remove contigs with more individuals exhibiting low coverage or genotyping errors [10]
        19 rad_haplotyper       -h      50              #Remove contigs with greater NumHaplotypes-NumSNPs. [100]
        19 rad_haplotyper       -z      0.2             #Remove up to this proportion or number of reads when testing for paralogs.  The more real variation in your data set, the greater this number will be. (<1) or number (>=1) of reads. [0.1]
        19 rad_haplotyper       -m      0.25            #Keep loci with a greater proportion of haplotyped individuals  [0.5]
```

5. Set up `fltrVCF.sbatch` to work with your system

6. Run `fltrVCF.sbatch`
```bash
sbatch fltrVCF.sbatch
```

7. Evaluate `slurm*out` and `*.pdf` that are created, adjust settings as required and rerun `fltrVCF.sbatch` . Save rad_haplotyper for later unless everything looks fantastic.

8. When `fltrVCF.sbatch` is done, run `fltrVCFstats2.sbatch` then evaluate `slurm*out` and `*pdf`.  Adjust settings and rerun 'fltrVCF.sbatch'

```bash
#the two arguments passed to fltrVCFstats2 are number of threads and the file prefix on all of the vcf files to evaluate
sbatch fltrVCFstats2.sbatch 20 PIRE_SiganusSpinus.K.10.10
```

9. Evaluate radhaplotyper output.  We will use this only to remove contigs from consideration - primarily potential paralogs.

10. Final goal is to make a list of contigs that were well represented and highly likely to be viable.  Those that made it through filter 041 minus those identified as paralogous by rad_haplotyper.


## Evaluating Output FILTER 30, Positions to be Considered

In the `slurm*out`, you should see some R feedback like this:

```bash
Wed Oct 9 23:30:49 CDT 2019 ---------------------------FILTER30: Remove Positions From Consideration -----------------------------
── Attaching packages ─────────────────────────────────────── tidyverse 1.2.1 ──
✔ ggplot2 3.2.0     ✔ purrr   0.3.2
✔ tibble  2.1.3     ✔ dplyr   0.8.3
✔ tidyr   1.0.0     ✔ stringr 1.4.0
✔ readr   1.3.1     ✔ forcats 0.4.0
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()

Attaching package: ‘gridExtra’

The following object is masked from ‘package:dplyr’:

    combine

There were 20 warnings (use warnings() to see them)
Warning messages:
1: Removed 6326759 rows containing missing values (stat_boxplot).
2: Removed 1 rows containing missing values (geom_segment).
3: Removed 1 rows containing missing values (geom_segment).
4: Removed 7049807 rows containing missing values (stat_boxplot).
5: Removed 1 rows containing missing values (geom_segment).
6: Removed 1 rows containing missing values (geom_segment).
Warning messages:
1: Removed 102548 rows containing non-finite values (stat_bin).
2: Removed 2 rows containing missing values (geom_bar).
Warning messages:
1: Removed 6932 rows containing non-finite values (stat_bin).
2: Removed 2 rows containing missing values (geom_bar).
Warning messages:
1: Removed 5889 rows containing non-finite values (stat_bin).
2: Removed 2 rows containing missing values (geom_bar).
null device
          1
 Plots output to PIRE_SiganusSpinus.Z.15.15.Fltr30.1.prefltr.plots.pdf
── Attaching packages ─────────────────────────────────────── tidyverse 1.2.1 ──
✔ ggplot2 3.2.0     ✔ purrr   0.3.2
✔ tibble  2.1.3     ✔ dplyr   0.8.3
✔ tidyr   1.0.0     ✔ stringr 1.4.0
✔ readr   1.3.1     ✔ forcats 0.4.0
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()

Attaching package: ‘gridExtra’

The following object is masked from ‘package:dplyr’:

    combine

Warning messages:
1: Removed 1180725 rows containing missing values (stat_boxplot).
2: Removed 1 rows containing missing values (geom_segment).
3: Removed 1 rows containing missing values (geom_segment).
4: Removed 1209807 rows containing missing values (stat_boxplot).
5: Removed 1 rows containing missing values (geom_segment).
6: Removed 1 rows containing missing values (geom_segment).
Warning messages:
1: Removed 14495 rows containing non-finite values (stat_bin).
2: Removed 2 rows containing missing values (geom_bar).
Warning messages:
1: Removed 591 rows containing non-finite values (stat_bin).
2: Removed 2 rows containing missing values (geom_bar).
Warning messages:
1: Removed 767 rows containing non-finite values (stat_bin).
2: Removed 2 rows containing missing values (geom_bar).
null device
          1
 Plots output to PIRE_SiganusSpinus.Z.15.15.Fltr30.1.postfltr.plots.pdf
```

The pre and post filter `*pdf` are what you're interested in.

Pre filtering, the mean depth of coverage per site for each contig is plotted against position to help you determine where to set the position cutoffs.  Below, things start out poorly at position 1 but quickly improve by position 9, then go south just before position 150, before recovering around pos 250.  I targeted 9-141.
![image](https://user-images.githubusercontent.com/12803659/66541702-d1d82000-eaf5-11e9-82ef-5fd915753532.png)

There are additional plots depicting SNPs, insertions, and deletions by position also to help determine good cutoffs.

Post filter:
![image](https://user-images.githubusercontent.com/12803659/66541972-c9ccb000-eaf6-11e9-92f6-3163e788f1ad.png)
![image](https://user-images.githubusercontent.com/12803659/66542002-dd781680-eaf6-11e9-8303-89b014f95b36.png)
![image](https://user-images.githubusercontent.com/12803659/66542026-f254aa00-eaf6-11e9-803b-324a41957312.png)
![image](https://user-images.githubusercontent.com/12803659/66542036-00a2c600-eaf7-11e9-88ef-0128240b6562.png)

## Evaluating Output FILTER 14, Min Depth for Genotype to be Considered

Most of the output you are interested in is in `slurm*out`.  I started out conservative with a setting of 5, but you might want to bump it up to 10.

```bash
                                                   Histogram of genotype depth before FILTER14

                       4e+07 +--------------------------------------------------------------------------------------+
                             |****            +                 +                +                 +                |
                             |   *                     'genotypedepthbefore' using (bin($1,binwidth)):(1.0) ******* |
                     3.5e+07 |-+ *                                                                                +-|
                             |   *                                                                                  |
                             |   *                                                                                  |
                       3e+07 |-+ *                                                                                +-|
                             |   *                                                                                  |
                             |   *                                                                                  |
                     2.5e+07 |-+ *                                                                                +-|
                             |   *                                                                                  |
                       2e+07 |-+ *                                                                                +-|
                             |   *                                                                                  |
                             |   *                                                                                  |
                     1.5e+07 |-+ *                                                                                +-|
                             |   *                                                                                  |
                             |   *                                                                                  |
                       1e+07 |-+ *                                                                                +-|
                             |   *                                                                                  |
                             |   *                                                                                  |
                       5e+06 |-+ *                                                                                +-|
                             |   *                                                                                  |
                             |   ***************************************************************************        |
                           0 +--------------------------------------------------------------------------------------+
                             0                20                40               60                80              100
                                                                      Depth

^L
                                     Scatter plot of mean depth per site before FILTER14.

       100 +--------------------------------------------------------------------------------------------------------+
           |            +        **  +            +             +            +            +            +            |
           |                      **                                                  'genotypedepthbefore'    *    |
           |                       **                                                                               |
           |                        **                                                                              |
        80 |-+                       **                                                                           +-|
           |                          **                                                                            |
           |                           **                                                                           |
           |                            **                                                                          |
        60 |-+                           **                                                                       +-|
           |                              **                                                                        |
           |                               ***                                                                      |
           |                                 **                                                                     |
           |                                  **                                                                    |
        40 |-+                                 ***                                                                +-|
           |                                     **                                                                 |
           |                                      **                                                                |
           |                                       ***                                                              |
        20 |-+                                       **                                                           +-|
           |                                          ***                                                           |
           |                                            ***                                                         |
           |                                              ***                                                       |
           |            +            +            +         *****            +            +            +            |
         0 +--------------------------------------------------------------------------------------------------------+
           0          1e+07        2e+07        3e+07         4e+07        5e+07        6e+07        7e+07        8e+07
                                                           Genotype

     vcftools --vcf PIRE_SiganusSpinus.Z.15.15.Fltr30.1.vcf --minDP 5 --recode --recode-INFO-all --out PIRE_SiganusSpinus.Z.15.15.Fltr14.2.recode.vcf 2> /dev/null
        Sites remaining:        1934863
        Contigs remaining:      14916
```
## Evaluating Output FILTER 16, Individuals with Too Much Missing Data


```bash
                         Histogram of % missing data per individual before filter. Bars to the left are desireable.

                       8 +------------------------------------------------------------------------------------------+
                         |                 +                 +                  +                 +               **|
                         |                                    'totalmissing' using (bin($1,binwidth)):(1.0) ********|
                       7 |-+                                                                                      **|
                         |                                                                                        **|
                         |                                                                                        **|
                       6 |-+                                                                                      **|
                         |                                                                                        **|
                         |                                                                                        **|
                       5 |-+                                                                                      **|
                         |                                                                                        **|
                       4 |-+ **                                                                                   **|
                         |   **                                                                                   **|
                         |   **                                                                                   **|
                       3 |-+****                                                                                  **|
                         |  ****                                                                                  **|
                         |  ****                                                                                  **|
                       2 |-*****        ***********                                                               **|
                         | *****        *         *                                                               **|
                         | *****        *         *                                                               **|
                       1 |-**************         ******************************************************************|
                         | *****   *    *         *           *            *               *         *      *  * ***|
                         | *****   *    *  +      *          +*            *    +          *      +  *      *  * ***|
                       0 +------------------------------------------------------------------------------------------+
                         0                0.2               0.4                0.6               0.8                1
                                                             % missing genotypes

^L
                                              Scatter plot of % missing data per individual.

                     1 +--------------------------------------------------------------------------------------------+
                       |          +           +        * * *        +           +          +          +           + |
                       |                                      *                                 'imiss.dat'    *    |
                       |                                                                                            |
                       |                                                                                            |
                   0.8 |-+                                      *                                                 +-|
                       |                                                                                            |
                       |                                                                                            |
                       |                                          *                                                 |
                   0.6 |-+                                                                                        +-|
                       |                                                                                            |
                       |                                                                                            |
                       |                                                                                            |
                       |                                            *                                               |
                   0.4 |-+                                             *                                          +-|
                       |                                                                                            |
                       |                                                                                            |
                       |                                                                                            |
                   0.2 |-+                                                                                        +-|
                       |                                                 * *                                        |
                       |                                                     *                                      |
                       |                                                        * * *                               |
                       |          +           +          +          +           +      * * * *  * * * *  * * * *  + |
                     0 +--------------------------------------------------------------------------------------------+
                       0          5           10         15         20          25         30         35          40
                                                                Individual

 Individuals with too much missing data:
INDV
PIRE2019-Ssp-A-Atu_010-Plate1Pool5Seq1-1E-L4
PIRE2019-Ssp-A-Atu_015-Plate1Pool8Seq1-1H-L4
PIRE2019-Ssp-A-Atu_017-Plate1Pool6Seq1-1F-L4
PIRE2019-Ssp-A-Atu_021-Plate1Pool3Seq1-1C-L4
PIRE2019-Ssp-A-Atu_024-Plate1Pool7Seq1-1G-L4
PIRE2019-Ssp-A-Atu_028-Plate1Pool2Seq1-1B-L4
PIRE2019-Ssp-A-Atu_034-Plate1Pool4Seq1-1D-L4
PIRE2019-Ssp-A-Atu_037-Plate1Pool1Seq1-1A-L4
PIRE2019-Ssp-A-Atu_039-Plate1Pool1Seq1-1A-L4
PIRE2019-Ssp-A-Atu_041-Plate1Pool8Seq1-1H-L4
PIRE2019-Ssp-A-Atu_043-Plate1Pool6Seq1-1F-L4
PIRE2019-Ssp-A-Atu_044-Plate1Pool7Seq1-1G-L4
PIRE2019-Ssp-A-Atu_045-Plate1Pool2Seq1-1B-L4
PIRE2019-Ssp-A-Atu_046-Plate1Pool5Seq1-1E-L4
PIRE2019-Ssp-A-Atu_048-Plate1Pool4Seq1-1D-L4
PIRE2019-Ssp-A-Atu_050-Plate1Pool3Seq1-1C-L4
PIRE2019-Ssp-C-Gub_037-Plate1Pool4Seq1-2E-L4
PIRE2019-Ssp-C-Gub_038-Plate1Pool3Seq1-2D-L4
PIRE2019-Ssp-C-Gub_052-Plate1Pool7Seq1-2H-L4
     vcftools --vcf PIRE_SiganusSpinus.Z.15.15.Fltr14.2.recode.vcf --remove PIRE_SiganusSpinus.Z.15.15.Fltr16.3.lowDP-2.indv --recode --recode-INFO-all --out PIRE_SiganusSpinus.Z.15.15.Fltr16.3.recode.vcf 2> /dev/null
        Sites remaining:        1934863
        Contigs remaining:      14916

sed: couldn't close ./sedckyku8: Permission denied

```

I don't know what's going on with that sed error.  It's relatively new... but doesn't seem to be causing a failure of the filter

## Evaluating Output FILTER 31,  Contig Length

This is another new filter where I used R instead of gnu plot, so the `*pdf` will be of interest. The `slurm*out` should look something like this:

```bash
── Attaching packages ─────────────────────────────────────── tidyverse 1.2.1 ──
✔ ggplot2 3.2.0     ✔ purrr   0.3.2
✔ tibble  2.1.3     ✔ dplyr   0.8.3
✔ tidyr   1.0.0     ✔ stringr 1.4.0
✔ readr   1.3.1     ✔ forcats 0.4.0
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()

Attaching package: ‘gridExtra’

The following object is masked from ‘package:dplyr’:

    combine

`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
`geom_smooth()` using method = 'gam' and formula 'y ~ s(x, bs = "cs")'
null device
          1
 Plots output to PIRE_SiganusSpinus.Z.15.15.Fltr31.4.ldepth.mean.contigs.pdf
```

This pdf starts you with the filtered result:
![image](https://user-images.githubusercontent.com/12803659/66542292-d43b7980-eaf7-11e9-986d-47d764ff319d.png)

Then it shows you different views of the data before applying the filter to help you identify better settings:
![image](https://user-images.githubusercontent.com/12803659/66542327-eddcc100-eaf7-11e9-9fd7-286a4dcb9e18.png)

## Evaluating Output FILTER 041, Contig Coverage

This filter is for the mean of mean depths by site for each contig.  I removed contigs  <=1st percentile and >=99th percentile.

This is another new filter where I used R instead of gnu plot, so the `*pdf` will be of interest. The `slurm*out` should look something like this:

```bash
── Attaching packages ─────────────────────────────────────── tidyverse 1.2.1 ──
✔ ggplot2 3.2.0     ✔ purrr   0.3.2
✔ tibble  2.1.3     ✔ dplyr   0.8.3
✔ tidyr   1.0.0     ✔ stringr 1.4.0
✔ readr   1.3.1     ✔ forcats 0.4.0
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()

Attaching package: ‘gridExtra’

The following object is masked from ‘package:dplyr’:

    combine

Warning messages:
1: Removed 2 rows containing missing values (geom_bar).
2: Removed 2 rows containing missing values (geom_bar).
null device
          1
 Plots output to PIRE_SiganusSpinus.Z.15.15.Fltr041.5.ldepth.mean.contigs.pdf
```

The pdf here is pretty straight foward, 1 page of output, first histogram is pre filter and second is post filter:
![image](https://user-images.githubusercontent.com/12803659/66542466-43b16900-eaf8-11e9-9bc9-3ecce15a91aa.png)

## Evaluating Output FILTER 05, Sites With Too Much Missing Data

The purpose of this filter is to try to clean things up a bit for rad_haplotyper, although, it didn't do much of that for me here.

The relevant output is in `slurm*out`:

```bash
                                          Histogram of the proportion of genotypes that are missing

                 1.4e+06 +------------------------------------------------------------------------------------------+
                         |                 +                 +                  +                 +                 |
                         |                                 'sitemissingness' using (bin($1,binwidth)):(1.0) ******* |
                 1.2e+06 |*****                                                                                   +-|
                         |    *                                                                                     |
                         |    *                                                                                     |
                         |    *                                                                                     |
                   1e+06 |-+  *                                                                                   +-|
                         |    *                                                                                     |
                         |    *                                                                                     |
                  800000 |-+  *                                                                                   +-|
                         |    *                                                                                     |
                         |    *                                                                                     |
                  600000 |-+  *                                                                                   +-|
                         |    *                                                                                     |
                         |    *                                                                                     |
                  400000 |-+  *                                                                                   +-|
                         |    *****                                                                                 |
                         |    *   *                                                                                 |
                         |    *   *                                                                                 |
                  200000 |-+  *   *                                                                               +-|
                         |    *   ******                                                                            |
                         |    *   *    *****                 +                  +                 +                 |
                       0 +------------------------------------------------------------------------------------------+
                         0                0.2               0.4                0.6               0.8                1
                                                             Proportion Missing

                                             Scatter plot of the proportion of genotypes that are missing.

                                1 +---------------------------------------------------------------------------------+
                                  |       +       +        +       +       +       +       +        +       +       |
                                  |                                                       'sitemissingness'    *    |
                                  |                                                                                 |
                                  |                                                                                 |
                              0.8 |-+                                                                             +-|
                                  |                                                                                 |
                                  |                                                                                 |
                                  |*                                                                                |
                              0.6 |*+                                                                             +-|
                                  |*                                                                                |
                                  |*                                                                                |
                                  |*                                                                                |
                                  |*                                                                                |
                              0.4 |**                                                                             +-|
                                  | *                                                                               |
                                  | *                                                                               |
                                  |                                                                                 |
                              0.2 |-**                                                                            +-|
                                  |  ***                                                                            |
                                  |    ******                                                                       |
                                  |         ****************                                                        |
                                  |       +       +        ***************************     +        +       +       |
                                0 +---------------------------------------------------------------------------------+
                                  0     200000  400000   600000  800000  1e+06  1.2e+06 1.4e+06  1.6e+06 1.8e+06  2e+06
                                                                         Site

     vcftools --vcf PIRE_SiganusSpinus.Z.15.15.Fltr041.5.vcf --max-missing 0.5 --recode --recode-INFO-all --out PIRE_SiganusSpinus.Z.15.15.Fltr05.1.recode.vcf 2> /dev/null
        Sites remaining:        1781900
        Contigs remaining:      13490
```

## Evaluating FILTER 07, Remove Monomorphic Sites

The only purpose of this filter is to remove monomorphic sites so that rad_haplotyper will finish this century. With it set at 0, it will remove sites with no called alternate alleles.

The `slurm*out` has the relevant output
```bash
vcffilter -s -f "AC > 0 & AN - AC > 0" PIRE_SiganusSpinus.Z.15.15.Fltr05.1.recode.vcf > PIRE_SiganusSpinus.Z.15.15.Fltr07.2.vcf
        Sites remaining:        17889
        Contigs remaining:      9071
```
## Evaluating FILTER 19, rad_haplotyper
The goal here is to develop a list of contigs to remove from the Filter32 vcf and then use that list of PASSING contigs to filter the reference genome to obtain our sequences for probe development.

---
### Contamination
You should review the `PIRE_SiganusSpinus.Z.15.15.Fltr19.ind_stats.out` file. You should be looking for individuals that are outliers.  This is a way to identify contamination, e.g., and indiviudal with an excessive amount of paralogs.

```bash
less -S PIRE_SiganusSpinus.Z.15.15.Fltr19.ind_stats.out 

Ind     Poss_Paralogs   Low_Coverage/Errors     Miss_Genotype   Total_Failed    Total_Loci      Prop_Success
PIRE2019-Ssp-C-Gub_002-Plate1Pool3Seq1-2D-L4.15.15      35      15      1012    1062    9071    0.882923602689891
PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.15.15     12      13      39      64      9071    0.992944548561349
PIRE2019-Ssp-C-Gub_006-Plate1Pool11Seq1-3D-L4.15.15     4       10      3       17      9071    0.998125895711608
PIRE2019-Ssp-C-Gub_012-Plate1Pool9Seq1-3B-L4.15.15      9       7       3       19      9071    0.997905412854151
PIRE2019-Ssp-C-Gub_019-Plate1Pool10Seq1-3C-L4.15.15     10      10      8       28      9071    0.99691323999559
PIRE2019-Ssp-C-Gub_025-Plate1Pool9Seq1-3B-L4.15.15      11      7       7       25      9071    0.997243964281777
PIRE2019-Ssp-C-Gub_028-Plate1Pool8Seq1-3A-L4.15.15      2       7       2949    2958    9071    0.673905853819865
PIRE2019-Ssp-C-Gub_033-Plate1Pool2Seq1-2C-L4.15.15      16      3       26      45      9071    0.995039135707199
PIRE2019-Ssp-C-Gub_034-Plate1Pool1Seq1-2B-L4.15.15      10      7       13      30      9071    0.996692757138133
PIRE2019-Ssp-C-Gub_036-Plate1Pool7Seq1-2H-L4.15.15      5       13      778     796     9071    0.912247822731783
PIRE2019-Ssp-C-Gub_045-Plate1Pool2Seq1-2C-L4.15.15      13      8       85      106     9071    0.988314408554735
PIRE2019-Ssp-C-Gub_063-Plate1Pool1Seq1-2B-L4.15.15      10      1       1       12      9071    0.998677102855253
PIRE2019-Ssp-C-Gub_068-Plate1Pool6Seq1-2G-L4.15.15      7       10      27      44      9071    0.995149377135928
PIRE2019-Ssp-C-Gub_071-Plate1Pool5Seq1-2F-L4.15.15      9       4       65      78      9071    0.991401168559145
PIRE2019-Ssp-C-Gub_075-Plate1Pool13Seq1-3F-L4.15.15     13      7       48      68      9071    0.992503582846434
PIRE2019-Ssp-C-Gub_077-Plate1Pool4Seq1-2E-L4.15.15      4       13      3425    3442    9071    0.62054900231507
PIRE2019-Ssp-C-Gub_082-Plate1Pool5Seq1-2F-L4.15.15      11      6       18      35      9071    0.996141549994488
PIRE2019-Ssp-C-Gub_085-Plate1Pool8Seq1-3A-L4.15.15      7       7       46      60      9071    0.993385514276265
PIRE2019-Ssp-C-Gub_086-Plate1Pool12Seq1-3E-L4.15.15     6       13      937     956     9071    0.894609194135156
PIRE2019-Ssp-C-Gub_087-Plate1Pool14Seq1-3G-L4.15.15     4       4       4       12      9071    0.998677102855253
PIRE2019-Ssp-C-Gub_096-Plate1Pool6Seq1-2G-L4.15.15      7       8       17      32      9071    0.996472274280675
```

I don't see evidence of contamination (large prop of contigs possible paralogs) and we've filtered for the other issues that rad_haplotyper identifies already.  If you do have an issue here, remove the individual from the TotalRawSNPs*vcf file and filter again.

---

### Paralog Exploration

We are primarily interested in the `PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out` file and will use this to remove some loci from consideration due to possible paralogs.

First, make sure your settings were correct.  Contigs should only be removed because of "Complex" variants or an extreme amount of missing data.
```bash
grep 'FILTERED' PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | less -S

dDocent_Contig_10029    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_10320    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_10330    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_1072     -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_10868    10      4       4       19      0.211   FILTERED        14      1       0       Missing data
dDocent_Contig_11042    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11054    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11076    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11122    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11219    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11243    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11247    -       -       -       -       -       FILTERED        0       0       0       Complex
dDocent_Contig_11334    -       -       -       -       -       FILTERED        0       0       0       Complex

grep 'FILTERED' PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | cut -f11 | sort | uniq -c
    585 Complex
     14 Missing data
```

If you don't get something like the above, adjust your thresholds and run again.

Second, evaluate possible paralogs.  The big question is where to draw the line.  0 possible paralogs, 1, 2... ?
```bash
tail -n+3 PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | cut -f8 | sort -n | uniq -c
  24871 0
   1686 1
    166 2
     23 3
     14 4
     14 5
     12 6
      4 7
      5 8
      2 9
      1 10
      4 11
      2 12
      1 14
      2 15

grep 'PASSED' PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | sort -nrk8 | less -S

dDocent_Contig_55040    11      3       7       19      0.368   PASSED  12      0       0
dDocent_Contig_14387    14      6       6       19      0.316   PASSED  12      0       1
dDocent_Contig_87409    11      9       7       19      0.368   PASSED  11      0       1
dDocent_Contig_7382     7       5       8       19      0.421   PASSED  11      0       0
dDocent_Contig_4780     10      7       8       19      0.421   PASSED  11      0       0
dDocent_Contig_23173    9       6       8       19      0.421   PASSED  11      0       0
dDocent_Contig_6904     6       6       9       19      0.474   PASSED  10      0       0
dDocent_Contig_61565    11      6       9       19      0.474   PASSED  9       0       1
dDocent_Contig_53074    7       6       9       19      0.474   PASSED  9       0       1
dDocent_Contig_94038    5       8       10      19      0.526   PASSED  8       0       1
dDocent_Contig_90092    4       4       9       19      0.474   PASSED  8       0       2
dDocent_Contig_74764    8       6       8       19      0.421   PASSED  8       0       3
dDocent_Contig_21249    7       7       11      19      0.579   PASSED  8       0       0
dDocent_Contig_10501    6       7       11      19      0.579   PASSED  8       0       0
dDocent_Contig_89914    10      7       9       19      0.474   PASSED  7       0       3
dDocent_Contig_59476    5       6       8       19      0.421   PASSED  7       0       4
dDocent_Contig_45198    7       6       11      19      0.579   PASSED  7       1       0
dDocent_Contig_19199    11      6       12      19      0.632   PASSED  7       0       0
dDocent_Contig_93726    6       2       11      19      0.579   PASSED  6       0       2
dDocent_Contig_93178    4       6       7       19      0.368   PASSED  6       1       5
dDocent_Contig_92245    10      7       8       19      0.421   PASSED  6       1       4
dDocent_Contig_90243    12      7       11      19      0.579   PASSED  6       0       2
dDocent_Contig_8985     11      8       11      19      0.579   PASSED  6       0       2
dDocent_Contig_85804    5       9       12      19      0.632   PASSED  6       0       1
dDocent_Contig_85577    20      10      11      19      0.579   PASSED  6       0       2
dDocent_Contig_85191    22      12      8       19      0.421   PASSED  6       0       5
dDocent_Contig_5364     4       3       12      19      0.632   PASSED  6       1       0
dDocent_Contig_23257    15      10      10      19      0.526   PASSED  6       1       2
dDocent_Contig_20864    5       4       13      19      0.684   PASSED  6       0       0
dDocent_Contig_1929     8       9       13      19      0.684   PASSED  6       0       0
dDocent_Contig_8940     9       4       12      19      0.632   PASSED  5       1       1
dDocent_Contig_85421    14      8       7       19      0.368   PASSED  5       4       3
dDocent_Contig_82427    6       7       13      19      0.684   PASSED  5       0       1
dDocent_Contig_78152    7       3       11      19      0.579   PASSED  5       0       3
dDocent_Contig_72947    11      8       12      19      0.632   PASSED  5       0       2
dDocent_Contig_55605    6       9       12      19      0.632   PASSED  5       0       2
dDocent_Contig_4923     12      9       10      19      0.526   PASSED  5       0       4
dDocent_Contig_42088    4       5       14      19      0.737   PASSED  5       0       0
dDocent_Contig_39137    3       4       14      19      0.737   PASSED  5       0       0
dDocent_Contig_3275     15      4       13      19      0.684   PASSED  5       1       0
dDocent_Contig_24205    7       5       12      19      0.632   PASSED  5       0       2
dDocent_Contig_23614    10      5       10      19      0.526   PASSED  5       3       1
dDocent_Contig_1936     9       8       14      19      0.737   PASSED  5       0       0
dDocent_Contig_1442     5       8       14      19      0.737   PASSED  5       0       0
dDocent_Contig_94359    24      15      11      19      0.579   PASSED  4       1       3
dDocent_Contig_9152     4       5       14      19      0.737   PASSED  4       0       1
dDocent_Contig_85849    3       5       11      19      0.579   PASSED  4       0       4
dDocent_Contig_82052    5       7       12      19      0.632   PASSED  4       0       3
dDocent_Contig_79423    6       6       13      19      0.684   PASSED  4       0       2
dDocent_Contig_71535    7       6       13      19      0.684   PASSED  4       1       1
dDocent_Contig_71105    7       9       13      19      0.684   PASSED  4       0       2
dDocent_Contig_64129    5       7       13      19      0.684   PASSED  4       0       2
dDocent_Contig_60335    3       5       15      19      0.789   PASSED  4       0       0
dDocent_Contig_51426    4       5       11      19      0.579   PASSED  4       0       4
dDocent_Contig_33983    9       9       14      19      0.737   PASSED  4       0       1
dDocent_Contig_3307     11      3       10      19      0.526   PASSED  4       1       4
dDocent_Contig_26272    10      7       13      19      0.684   PASSED  4       1       1
dDocent_Contig_21773    7       5       15      19      0.789   PASSED  4       0       0
dDocent_Contig_88477    16      14      9       19      0.474   PASSED  3       1       6
dDocent_Contig_86710    10      12      15      19      0.789   PASSED  3       0       1
dDocent_Contig_86392    8       12      13      19      0.684   PASSED  3       1       2
dDocent_Contig_82786    9       12      16      19      0.842   PASSED  3       0       0
dDocent_Contig_82422    6       6       11      19      0.579   PASSED  3       0       5
dDocent_Contig_82344    6       4       13      19      0.684   PASSED  3       0       3
dDocent_Contig_78314    15      8       11      19      0.579   PASSED  3       2       3
dDocent_Contig_77997    5       8       15      19      0.789   PASSED  3       0       1
dDocent_Contig_75372    7       8       13      19      0.684   PASSED  3       0       3
dDocent_Contig_73265    18      9       12      19      0.632   PASSED  3       2       2
dDocent_Contig_71794    8       9       14      19      0.737   PASSED  3       0       2
dDocent_Contig_40667    9       6       14      19      0.737   PASSED  3       1       1
dDocent_Contig_33835    10      12      16      19      0.842   PASSED  3       0       0
dDocent_Contig_31418    9       9       11      19      0.579   PASSED  3       0       5
dDocent_Contig_31391    6       7       11      19      0.579   PASSED  3       0       5
dDocent_Contig_31270    14      3       15      19      0.789   PASSED  3       0       1
dDocent_Contig_24423    13      9       12      19      0.632   PASSED  3       2       2
dDocent_Contig_21075    8       5       14      19      0.737   PASSED  3       0       2
dDocent_Contig_18786    5       7       16      19      0.842   PASSED  3       0       0
dDocent_Contig_152      11      9       16      19      0.842   PASSED  3       0       0
dDocent_Contig_13997    6       6       15      19      0.789   PASSED  3       0       1
dDocent_Contig_13029    2       3       14      19      0.737   PASSED  3       0       2
dDocent_Contig_12853    4       6       14      19      0.737   PASSED  3       0       2
dDocent_Contig_97564    8       12      15      19      0.789   PASSED  2       1       1
dDocent_Contig_97124    5       6       16      19      0.842   PASSED  2       0       1
dDocent_Contig_9599     5       6       17      19      0.895   PASSED  2       0       0
dDocent_Contig_95398    3       4       15      19      0.789   PASSED  2       0       2
dDocent_Contig_94860    5       8       16      19      0.842   PASSED  2       0       1
dDocent_Contig_94385    6       13      14      19      0.737   PASSED  2       0       3
dDocent_Contig_94140    3       6       14      19      0.737   PASSED  2       0       3
dDocent_Contig_93974    5       7       17      19      0.895   PASSED  2       0       0
dDocent_Contig_93313    6       9       17      19      0.895   PASSED  2       0       0
dDocent_Contig_93080    7       10      16      19      0.842   PASSED  2       0       1
dDocent_Contig_92551    5       10      16      19      0.842   PASSED  2       0       1
dDocent_Contig_92052    3       4       14      19      0.737   PASSED  2       1       2
dDocent_Contig_90464    7       7       15      19      0.789   PASSED  2       0       2
dDocent_Contig_89949    10      11      16      19      0.842   PASSED  2       0       1
dDocent_Contig_89600    6       13      16      19      0.842   PASSED  2       0       1
dDocent_Contig_88373    9       16      16      19      0.842   PASSED  2       0       1
dDocent_Contig_88296    7       10      15      19      0.789   PASSED  2       0       2
dDocent_Contig_88189    4       4       15      19      0.789   PASSED  2       0       2
dDocent_Contig_87590    4       6       16      19      0.842   PASSED  2       0       1
dDocent_Contig_87110    4       7       14      19      0.737   PASSED  2       0       3
dDocent_Contig_85457    5       6       17      19      0.895   PASSED  2       0       0
dDocent_Contig_85333    8       15      16      19      0.842   PASSED  2       0       1
dDocent_Contig_85276    5       6       15      19      0.789   PASSED  2       0       2
dDocent_Contig_84690    6       9       14      19      0.737   PASSED  2       0       3
dDocent_Contig_83962    6       11      16      19      0.842   PASSED  2       0       1
dDocent_Contig_83740    3       5       16      19      0.842   PASSED  2       0       1
dDocent_Contig_83619    2       3       17      19      0.895   PASSED  2       0       0
dDocent_Contig_83401    9       15      17      19      0.895   PASSED  2       0       0
dDocent_Contig_82985    8       9       15      19      0.789   PASSED  2       0       2
dDocent_Contig_82920    4       7       16      19      0.842   PASSED  2       0       1
dDocent_Contig_82709    5       6       16      19      0.842   PASSED  2       0       1
dDocent_Contig_82453    9       8       16      19      0.842   PASSED  2       0       1
dDocent_Contig_82040    7       13      14      19      0.737   PASSED  2       2       1
dDocent_Contig_81893    5       11      16      19      0.842   PASSED  2       0       1
```

The expected number of paralogous individuals varies with the alternate allele frequency in truly paralogous contigs.  More variation, more paralogs.  You can cross-reference against `PIRE_SiganusSpinus.Z.15.15.Fltr19.haplo_dump.out`

```bash
less -S PIRE_SiganusSpinus.Z.15.15.Fltr19.haplo_dump.out

#search for the contig w/ most possible paralogs
          'dDocent_Contig_9361' => {
                                     'PIRE2019-Ssp-C-Gub_071-Plate1Pool5Seq1-2F-L4.15.15' => [
                                                                                               'GAC',
                                                                                               'AAC'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_036-Plate1Pool7Seq1-2H-L4.15.15' => [
                                                                                               'GGT',
                                                                                               'GAC'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_005-Plate1Pool12Seq1-3E-L4.15.15' => [
                                                                                                'GGT',
                                                                                                'GAC'
                                                                                              ],
                                     'PIRE2019-Ssp-C-Gub_087-Plate1Pool14Seq1-3G-L4.15.15' => [
                                                                                                'GAC',
                                                                                                'AAC'
                                                                                              ],
                                     'PIRE2019-Ssp-C-Gub_002-Plate1Pool3Seq1-2D-L4.15.15' => [
                                                                                               'GAT',
                                                                                               'GGT'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_006-Plate1Pool11Seq1-3D-L4.15.15' => [
                                                                                                'GAC',
                                                                                                'AAC'
                                                                                              ],
                                     'PIRE2019-Ssp-C-Gub_012-Plate1Pool9Seq1-3B-L4.15.15' => [
                                                                                               'GGT',
                                                                                               'GAC'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_028-Plate1Pool8Seq1-3A-L4.15.15' => [
                                                                                               'GAC'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_045-Plate1Pool2Seq1-2C-L4.15.15' => [
                                                                                               'GGT',
                                                                                               'AAC'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_033-Plate1Pool2Seq1-2C-L4.15.15' => [
                                                                                               'GGT',
                                                                                               'GAC'
                                                                                             ],
                                     'PIRE2019-Ssp-C-Gub_082-Plate1Pool5Seq1-2F-L4.15.15' => [
                                                                                               'AAC',
                                                                                               'GAC'
                                                                                             ]
                                   },
```

The high proportion of heterozygotes is highly suspicious.  

We can also cross reference against `PIRE_SiganusSpinus.Z.15.15.Fltr19.hap_log.out`.  This file has the information rad_haplotyper evaluated for each individual and each contig.

```bash
less -S PIRE_SiganusSpinus.Z.15.15.Fltr19.hap_log.out

#search for contig with most possible paralogs
dDocent_Contig_9361: Observed Haps:
$VAR1 = [
          'GGT',
          'GGT',
          'GAC',
          'GAC',
          'GGT',
          'GGT',
          'GGT',
          'GAT',
          'GAT',
          'GAT'
        ];
dDocent_Contig_9361: Unique Observed Haps:
$VAR1 = [
          'GGT',
          'GAC',
          'GAT'
        ];
dDocent_Contig_9361: Problem- trying to fix...
dDocent_Contig_9361: haplotype count threshold:2
dDocent_Contig_9361: Corrected Unique Observed Haps:
$VAR1 = [
          'GAT',
          'GGT'
        ];
dDocent_Contig_9361: Problem fixed
```

Looking at the above, my settings were probably too lenient (-z 0.2). I updated them to 0.1. Notice that there appears to be 3 haps, but one is removed from consideration. For this locus, it was ok, but maybe there's a few other loci where there are more than 2 haplotypes in and individual but the -z parameter is hiding them. Be sure to view more than 1 individual.

### More Paralog Exploration

Overall, I like the `*hap_log.out` file so I made a script to pull more information from it.  `rad_haplotyper` tries to rescue paralogous individuals with the aforementioned `-z` parameter, as well as evaluating all reads if a contig fails, and then applies the -z parameter again if necessary.  It also automatically passes SNPs as far as I know. Run `prlgStats.sbatch` to view these stats more closely (it's in the `fltrVCF/scripts` dir).  I'm particularly interested in filtering contigs with a high number of "rescued" individuals and reduces the need to experiment as much with the `-z` parameter.

```bash
sbatch prlgStats.sbatch 20 PIRE_SiganusSpinus.L.5.5

#sort by number of rescued individuals per contig
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | sort -nrk6 | less -S

#sort by num paralogs
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | sort -nrk9 | less -S

#count contigs with 1 SNP
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | cut -f5 | grep -c '[1-9][0-9]*'

#count contigs with >1 SNP
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | cut -f4 | grep -c '[1-9][0-9]*'

```

We can see the effect of applying filters for paralogs  and / or number of individuals rescued, note you may have to change the $10 to $9  and $7 to $6 (column identifiers) because of changes I've made to the `prlgStats.sbatch` script
```bash
#set min number of acceptable paralogous individuals
THRESHOLD=1
#set min number of acceptable rescued paralogous individuals (by the -z option)
THRESHOLDb=1

#count the number of contigs total in radhap output
tail -n+3 PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | wc -l
26807

# count number of non-complex contigs total
tail -n+3 PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | grep -vc Complex 
26222

# the number of contigs in your `*paralog_test_results.tsv` should be the same as the non-complex contigs above
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | wc -l
26222

#count the number of contigs passing paralog filters, change thresholds as neccessary
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 <= X {print ;}' | awk -v Y=$THRESHOLDb '$7 <= Y {print ;}' | wc -l
24702

THRESHOLD=0
THRESHOLDb=0
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 <= X {print ;}' | awk -v Y=$THRESHOLDb '$7 <= Y {print ;}' | wc -l
20173

THRESHOLD=0
THRESHOLDb=100
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 <= X {print ;}' | awk -v Y=$THRESHOLDb '$7 <= Y {print ;}' | wc -l
23391

THRESHOLD=1
THRESHOLDb=2
tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 <= X {print ;}' | awk -v Y=$THRESHOLDb '$7 <= Y {print ;}' | wc -l
25381
```

### Paralog Filter

We will use `awk` to start building a list of contigs to filter by slightly modifying the previous code

```bash
# list contigs with >= X paralogs and Y rescued paralogs, then save them into a variable, note you may have to change the column identifiers ($10, $7) because of changes I've made to the `prlgStats.sbatch` script since creating my `*paralog_test_results.tsv`
THRESHOLD=1
THRESHOLDb=2

tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 > X {print $1;}'  | less -S

tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v Y=$THRESHOLDb '$7> X {print $1;}' | less -S

PrlgCntgFltr1=$(tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v X=$THRESHOLD '$10 > X {print $1;}' )
PrlgCntgFltr2=$(tail -n+2 PIRE_SiganusSpinus.L.5.5.paralog_test_results.tsv | awk -v Y=$THRESHOLDb '$7> Y {print $1;}' )

# Make a third variable with the list of haps with extreme amounts of missing data
PrlgCntgFltr3=$(tail -n+3 PIRE_SiganusSpinus.L.5.5.Fltr19.stats.out | grep 'Missing data$' | cut -f1 )

# Concatenate the contigs in the PrlgCntgFltr vars and remove duplicates
cat <(echo $PrlgCntgFltr1 | tr " " "\n") <(echo $PrlgCntgFltr2 | tr " " "\n") <(echo $PrlgCntgFltr3 | tr " " "\n") | sort | uniq | sed -e 's/^/\^/g' -e 's/$/\$/g' > PIRE_SiganusSpinus.L.5.5.Fltr19.remove.contigs

# Get a list of the contigs that passed Filter 32 and remove those that failed Filter 19 and convert the contig list to a regex pattern
grep '^dDocent_Contig' PIRE_SiganusSpinus.L.5.5.Fltr32.6.vcf | cut -f1 | uniq | grep -vf PIRE_SiganusSpinus.L.5.5.Fltr19.remove.contigs | sed -e 's/^/\^>/g' -e 's/$/\$/g' > PIRE_SiganusSpinus.L.5.5.Fltr19.keep.contigs

```

###  Filter Reference to Create Fasta With Contigs for Probe Development

```bash
# Filter the reference genome for the contigs that passed filters, these are the contigs that will be considered for probe development
grep -A1 -f PIRE_SiganusSpinus.L.5.5.Fltr19.keep.contigs ../mkVCF2/reference.5.5.fasta | grep -v '^--$' > PIRE_SiganusSpinus.L.5.5.probes4development.fasta

# Remove loci with 24 single nuc repeats, 12 dinuc repeats, 8 trinuc repeats, or 6 tetranuc repeats
# See code block below for contents of microsat.motifs
cat PIRE_SiganusSpinus.L.5.5.probes4development.fasta | paste - - > PIRE_SiganusSpinus.L.5.5.probes4development.tsv

grep -vF -f microsat.motifs PIRE_SiganusSpinus.L.5.5.probes4development.tsv | tr "\t" "\n" > PIRE_SiganusSpinus.L.5.5.probes4development.noMSATS.fasta

```

`microsat.motifs`
```bash
GGGGGGGGGGGGGGGGGGGGGGGG
GAGAGAGAGAGAGAGAGAGAGAGA
GTGTGTGTGTGTGTGTGTGTGTGT
GCGCGCGCGCGCGCGCGCGCGCGC
AGAGAGAGAGAGAGAGAGAGAGAG
AAAAAAAAAAAAAAAAAAAAAAAA
ATATATATATATATATATATATAT
ACACACACACACACACACACACAC
TGTGTGTGTGTGTGTGTGTGTGTG
TATATATATATATATATATATATA
TTTTTTTTTTTTTTTTTTTTTTTT
TCTCTCTCTCTCTCTCTCTCTCTC
CGCGCGCGCGCGCGCGCGCGCGCG
CACACACACACACACACACACACA
CTCTCTCTCTCTCTCTCTCTCTCT
CCCCCCCCCCCCCCCCCCCCCCCC
GAAGAAGAAGAAGAAGAAGAAGAA
GTTGTTGTTGTTGTTGTTGTTGTT
GCCGCCGCCGCCGCCGCCGCCGCC
AGGAGGAGGAGGAGGAGGAGGAGG
ATTATTATTATTATTATTATTATT
ACCACCACCACCACCACCACCACC
TGGTGGTGGTGGTGGTGGTGGTGG
TAATAATAATAATAATAATAATAA
TCCTCCTCCTCCTCCTCCTCCTCC
CGGCGGCGGCGGCGGCGGCGGCGG
CAACAACAACAACAACAACAACAA
CTTCTTCTTCTTCTTCTTCTTCTT
GGAGGAGGAGGAGGAGGAGGAGGA
GATGATGATGATGATGATGATGAT
GTCGTCGTCGTCGTCGTCGTCGTC
GCGGCGGCGGCGGCGGCGGCGGCG
AGAAGAAGAAGAAGAAGAAGAAGA
AATAATAATAATAATAATAATAAT
ATCATCATCATCATCATCATCATC
ACGACGACGACGACGACGACGACG
TGATGATGATGATGATGATGATGA
TATTATTATTATTATTATTATTAT
TTCTTCTTCTTCTTCTTCTTCTTC
TCGTCGTCGTCGTCGTCGTCGTCG
CGACGACGACGACGACGACGACGA
CATCATCATCATCATCATCATCAT
CTCCTCCTCCTCCTCCTCCTCCTC
CCGCCGCCGCCGCCGCCGCCGCCG
GGTGGTGGTGGTGGTGGTGGTGGT
GACGACGACGACGACGACGACGAC
GTGGTGGTGGTGGTGGTGGTGGTG
GCAGCAGCAGCAGCAGCAGCAGCA
AGTAGTAGTAGTAGTAGTAGTAGT
AACAACAACAACAACAACAACAAC
ATGATGATGATGATGATGATGATG
ACAACAACAACAACAACAACAACA
TGTTGTTGTTGTTGTTGTTGTTGT
TACTACTACTACTACTACTACTAC
TTGTTGTTGTTGTTGTTGTTGTTG
TCATCATCATCATCATCATCATCA
CGTCGTCGTCGTCGTCGTCGTCGT
CACCACCACCACCACCACCACCAC
CTGCTGCTGCTGCTGCTGCTGCTG
CCACCACCACCACCACCACCACCA
GGCGGCGGCGGCGGCGGCGGCGGC
GAGGAGGAGGAGGAGGAGGAGGAG
GTAGTAGTAGTAGTAGTAGTAGTA
GCTGCTGCTGCTGCTGCTGCTGCT
AGCAGCAGCAGCAGCAGCAGCAGC
AAGAAGAAGAAGAAGAAGAAGAAG
ATAATAATAATAATAATAATAATA
ACTACTACTACTACTACTACTACT
TGCTGCTGCTGCTGCTGCTGCTGC
TAGTAGTAGTAGTAGTAGTAGTAG
TTATTATTATTATTATTATTATTA
TCTTCTTCTTCTTCTTCTTCTTCT
CGCCGCCGCCGCCGCCGCCGCCGC
CAGCAGCAGCAGCAGCAGCAGCAG
CTACTACTACTACTACTACTACTA
CCTCCTCCTCCTCCTCCTCCTCCT
GAAGGAAGGAAGGAAGGAAGGAAG
GTTGGTTGGTTGGTTGGTTGGTTG
GCCGGCCGGCCGGCCGGCCGGCCG
AGGGAGGGAGGGAGGGAGGGAGGG
AAAGAAAGAAAGAAAGAAAGAAAG
ATTGATTGATTGATTGATTGATTG
ACCGACCGACCGACCGACCGACCG
TGGGTGGGTGGGTGGGTGGGTGGG
TAAGTAAGTAAGTAAGTAAGTAAG
TTTGTTTGTTTGTTTGTTTGTTTG
TCCGTCCGTCCGTCCGTCCGTCCG
CGGGCGGGCGGGCGGGCGGGCGGG
CAAGCAAGCAAGCAAGCAAGCAAG
CTTGCTTGCTTGCTTGCTTGCTTG
CCCGCCCGCCCGCCCGCCCGCCCG
GGAGGGAGGGAGGGAGGGAGGGAG
GATGGATGGATGGATGGATGGATG
GTCGGTCGGTCGGTCGGTCGGTCG
GCGGGCGGGCGGGCGGGCGGGCGG
AATGAATGAATGAATGAATGAATG
ATCGATCGATCGATCGATCGATCG
ACGGACGGACGGACGGACGGACGG
TGAGTGAGTGAGTGAGTGAGTGAG
TATGTATGTATGTATGTATGTATG
TTCGTTCGTTCGTTCGTTCGTTCG
TCGGTCGGTCGGTCGGTCGGTCGG
CGAGCGAGCGAGCGAGCGAGCGAG
CATGCATGCATGCATGCATGCATG
CTCGCTCGCTCGCTCGCTCGCTCG
CCGGCCGGCCGGCCGGCCGGCCGG
GGTGGGTGGGTGGGTGGGTGGGTG
GACGGACGGACGGACGGACGGACG
GTGGGTGGGTGGGTGGGTGGGTGG
GCAGGCAGGCAGGCAGGCAGGCAG
AGTGAGTGAGTGAGTGAGTGAGTG
AACGAACGAACGAACGAACGAACG
ATGGATGGATGGATGGATGGATGG
ACAGACAGACAGACAGACAGACAG
TACGTACGTACGTACGTACGTACG
TTGGTTGGTTGGTTGGTTGGTTGG
TCAGTCAGTCAGTCAGTCAGTCAG
CGTGCGTGCGTGCGTGCGTGCGTG
CACGCACGCACGCACGCACGCACG
CTGGCTGGCTGGCTGGCTGGCTGG
CCAGCCAGCCAGCCAGCCAGCCAG
GGCGGGCGGGCGGGCGGGCGGGCG
GAGGGAGGGAGGGAGGGAGGGAGG
GTAGGTAGGTAGGTAGGTAGGTAG
GCTGGCTGGCTGGCTGGCTGGCTG
AGCGAGCGAGCGAGCGAGCGAGCG
AAGGAAGGAAGGAAGGAAGGAAGG
ATAGATAGATAGATAGATAGATAG
ACTGACTGACTGACTGACTGACTG
TGCGTGCGTGCGTGCGTGCGTGCG
TAGGTAGGTAGGTAGGTAGGTAGG
TTAGTTAGTTAGTTAGTTAGTTAG
TCTGTCTGTCTGTCTGTCTGTCTG
CAGGCAGGCAGGCAGGCAGGCAGG
CTAGCTAGCTAGCTAGCTAGCTAG
CCTGCCTGCCTGCCTGCCTGCCTG
GGGAGGGAGGGAGGGAGGGAGGGA
GAAAGAAAGAAAGAAAGAAAGAAA
GTTAGTTAGTTAGTTAGTTAGTTA
GCCAGCCAGCCAGCCAGCCAGCCA
AGGAAGGAAGGAAGGAAGGAAGGA
ATTAATTAATTAATTAATTAATTA
ACCAACCAACCAACCAACCAACCA
TGGATGGATGGATGGATGGATGGA
TAAATAAATAAATAAATAAATAAA
TTTATTTATTTATTTATTTATTTA
TCCATCCATCCATCCATCCATCCA
CGGACGGACGGACGGACGGACGGA
CAAACAAACAAACAAACAAACAAA
CTTACTTACTTACTTACTTACTTA
CCCACCCACCCACCCACCCACCCA
GGAAGGAAGGAAGGAAGGAAGGAA
GATAGATAGATAGATAGATAGATA
GTCAGTCAGTCAGTCAGTCAGTCA
GCGAGCGAGCGAGCGAGCGAGCGA
AGAAAGAAAGAAAGAAAGAAAGAA
AATAAATAAATAAATAAATAAATA
ATCAATCAATCAATCAATCAATCA
ACGAACGAACGAACGAACGAACGA
TGAATGAATGAATGAATGAATGAA
TTCATTCATTCATTCATTCATTCA
TCGATCGATCGATCGATCGATCGA
CGAACGAACGAACGAACGAACGAA
CATACATACATACATACATACATA
CTCACTCACTCACTCACTCACTCA
CCGACCGACCGACCGACCGACCGA
GGTAGGTAGGTAGGTAGGTAGGTA
GACAGACAGACAGACAGACAGACA
GTGAGTGAGTGAGTGAGTGAGTGA
GCAAGCAAGCAAGCAAGCAAGCAA
AGTAAGTAAGTAAGTAAGTAAGTA
AACAAACAAACAAACAAACAAACA
ATGAATGAATGAATGAATGAATGA
ACAAACAAACAAACAAACAAACAA
TGTATGTATGTATGTATGTATGTA
TACATACATACATACATACATACA
TTGATTGATTGATTGATTGATTGA
TCAATCAATCAATCAATCAATCAA
CGTACGTACGTACGTACGTACGTA
CTGACTGACTGACTGACTGACTGA
CCAACCAACCAACCAACCAACCAA
GGCAGGCAGGCAGGCAGGCAGGCA
GTAAGTAAGTAAGTAAGTAAGTAA
GCTAGCTAGCTAGCTAGCTAGCTA
AGCAAGCAAGCAAGCAAGCAAGCA
AAGAAAGAAAGAAAGAAAGAAAGA
ATAAATAAATAAATAAATAAATAA
ACTAACTAACTAACTAACTAACTA
TGCATGCATGCATGCATGCATGCA
TAGATAGATAGATAGATAGATAGA
TTAATTAATTAATTAATTAATTAA
TCTATCTATCTATCTATCTATCTA
CGCACGCACGCACGCACGCACGCA
CAGACAGACAGACAGACAGACAGA
CTAACTAACTAACTAACTAACTAA
CCTACCTACCTACCTACCTACCTA
GGGTGGGTGGGTGGGTGGGTGGGT
GAATGAATGAATGAATGAATGAAT
GTTTGTTTGTTTGTTTGTTTGTTT
GCCTGCCTGCCTGCCTGCCTGCCT
AGGTAGGTAGGTAGGTAGGTAGGT
AAATAAATAAATAAATAAATAAAT
ATTTATTTATTTATTTATTTATTT
ACCTACCTACCTACCTACCTACCT
TGGTTGGTTGGTTGGTTGGTTGGT
TAATTAATTAATTAATTAATTAAT
TCCTTCCTTCCTTCCTTCCTTCCT
CGGTCGGTCGGTCGGTCGGTCGGT
CAATCAATCAATCAATCAATCAAT
CTTTCTTTCTTTCTTTCTTTCTTT
CCCTCCCTCCCTCCCTCCCTCCCT
GGATGGATGGATGGATGGATGGAT
GATTGATTGATTGATTGATTGATT
GTCTGTCTGTCTGTCTGTCTGTCT
GCGTGCGTGCGTGCGTGCGTGCGT
AGATAGATAGATAGATAGATAGAT
AATTAATTAATTAATTAATTAATT
ATCTATCTATCTATCTATCTATCT
ACGTACGTACGTACGTACGTACGT
TGATTGATTGATTGATTGATTGAT
TATTTATTTATTTATTTATTTATT
TTCTTTCTTTCTTTCTTTCTTTCT
TCGTTCGTTCGTTCGTTCGTTCGT
CGATCGATCGATCGATCGATCGAT
CATTCATTCATTCATTCATTCATT
CCGTCCGTCCGTCCGTCCGTCCGT
GGTTGGTTGGTTGGTTGGTTGGTT
GACTGACTGACTGACTGACTGACT
GCATGCATGCATGCATGCATGCAT
AGTTAGTTAGTTAGTTAGTTAGTT
AACTAACTAACTAACTAACTAACT
ATGTATGTATGTATGTATGTATGT
ACATACATACATACATACATACAT
TGTTTGTTTGTTTGTTTGTTTGTT
TACTTACTTACTTACTTACTTACT
TTGTTTGTTTGTTTGTTTGTTTGT
TCATTCATTCATTCATTCATTCAT
CGTTCGTTCGTTCGTTCGTTCGTT
CACTCACTCACTCACTCACTCACT
CTGTCTGTCTGTCTGTCTGTCTGT
CCATCCATCCATCCATCCATCCAT
GGCTGGCTGGCTGGCTGGCTGGCT
GAGTGAGTGAGTGAGTGAGTGAGT
GTATGTATGTATGTATGTATGTAT
GCTTGCTTGCTTGCTTGCTTGCTT
AGCTAGCTAGCTAGCTAGCTAGCT
AAGTAAGTAAGTAAGTAAGTAAGT
ACTTACTTACTTACTTACTTACTT
TGCTTGCTTGCTTGCTTGCTTGCT
TAGTTAGTTAGTTAGTTAGTTAGT
TTATTTATTTATTTATTTATTTAT
TCTTTCTTTCTTTCTTTCTTTCTT
CGCTCGCTCGCTCGCTCGCTCGCT
CAGTCAGTCAGTCAGTCAGTCAGT
CTATCTATCTATCTATCTATCTAT
CCTTCCTTCCTTCCTTCCTTCCTT
GGGCGGGCGGGCGGGCGGGCGGGC
GAACGAACGAACGAACGAACGAAC
GTTCGTTCGTTCGTTCGTTCGTTC
GCCCGCCCGCCCGCCCGCCCGCCC
AGGCAGGCAGGCAGGCAGGCAGGC
AAACAAACAAACAAACAAACAAAC
ATTCATTCATTCATTCATTCATTC
ACCCACCCACCCACCCACCCACCC
TGGCTGGCTGGCTGGCTGGCTGGC
TAACTAACTAACTAACTAACTAAC
TTTCTTTCTTTCTTTCTTTCTTTC
TCCCTCCCTCCCTCCCTCCCTCCC
CGGCCGGCCGGCCGGCCGGCCGGC
CAACCAACCAACCAACCAACCAAC
CTTCCTTCCTTCCTTCCTTCCTTC
GGACGGACGGACGGACGGACGGAC
GATCGATCGATCGATCGATCGATC
GTCCGTCCGTCCGTCCGTCCGTCC
AGACAGACAGACAGACAGACAGAC
AATCAATCAATCAATCAATCAATC
ATCCATCCATCCATCCATCCATCC
ACGCACGCACGCACGCACGCACGC
TGACTGACTGACTGACTGACTGAC
TATCTATCTATCTATCTATCTATC
TTCCTTCCTTCCTTCCTTCCTTCC
TCGCTCGCTCGCTCGCTCGCTCGC
CGACCGACCGACCGACCGACCGAC
CATCCATCCATCCATCCATCCATC
CTCCCTCCCTCCCTCCCTCCCTCC
CCGCCCGCCCGCCCGCCCGCCCGC
GGTCGGTCGGTCGGTCGGTCGGTC
GACCGACCGACCGACCGACCGACC
GTGCGTGCGTGCGTGCGTGCGTGC
GCACGCACGCACGCACGCACGCAC
AGTCAGTCAGTCAGTCAGTCAGTC
AACCAACCAACCAACCAACCAACC
ATGCATGCATGCATGCATGCATGC
TGTCTGTCTGTCTGTCTGTCTGTC
TACCTACCTACCTACCTACCTACC
TTGCTTGCTTGCTTGCTTGCTTGC
TCACTCACTCACTCACTCACTCAC
CGTCCGTCCGTCCGTCCGTCCGTC
CACCCACCCACCCACCCACCCACC
CTGCCTGCCTGCCTGCCTGCCTGC
CCACCCACCCACCCACCCACCCAC
GGCCGGCCGGCCGGCCGGCCGGCC
GAGCGAGCGAGCGAGCGAGCGAGC
GTACGTACGTACGTACGTACGTAC
GCTCGCTCGCTCGCTCGCTCGCTC
AGCCAGCCAGCCAGCCAGCCAGCC
AAGCAAGCAAGCAAGCAAGCAAGC
ATACATACATACATACATACATAC
ACTCACTCACTCACTCACTCACTC
TGCCTGCCTGCCTGCCTGCCTGCC
TAGCTAGCTAGCTAGCTAGCTAGC
TTACTTACTTACTTACTTACTTAC
CGCCCGCCCGCCCGCCCGCCCGCC
CAGCCAGCCAGCCAGCCAGCCAGC
CTACCTACCTACCTACCTACCTAC
CCTCCCTCCCTCCCTCCCTCCCTC
```

## Evaluating FILTER 32: Remove Contigs with Excess Heterozygosity

This filter is meant to remove paralogs that I don't think rad_haplotyper can identify.  If every SNP is heterozygous in every individual in a contig, it is highly likely to be due to two or more invariant loci being mapped to the same contig.  For biallelic SNPs, we expect a max heterozygosity of 0.5.  I padded this and set my cutoff at 0.9.  It might be too lenient, but probably shouldn't go lower than 0.8 due to variation in sampling which may result in higher values.

I've begun troubleshooting this, fixed.

pull the newest version of fltrVCF
fltrVCF.bash updated

# Development of Filter for CV of Coverage

Why CV? Because each contig has a different mean depth of coverage and CV is standardized by the mean, allowing for apples-to-apples comparisons that aren't biased by cvg.

I'm modifying Fltr041, which filters for extreme coverage. It is important to note that we are stuck with working position by position for these filters due to the nature of the VCF.  Consequently, we have the mean and var of coverage depth across individuals per position to work with but need to scale up to the contig.

## Filter041 for mean cvg by contig 
This filter works by taking the contig-mean of mean cvgs per position, then removing contigs that fall on the fringes of the frequency distribution.  When I check it by eye, the distribution is approximately normal and the extreme contigs have consistently low or high cvg across all positions, so this seems to work well.

## Filtering for the coefficient of variation in cvg by contig
Here, rather than calculating the mean of the mean cvg by position, we calculate the *_coefficient of variation of the mean of the mean cvg by position_*.  That's a mouthful.  When I check it by eye, the distribution pinned against zero (sort of like a Poisson dist).  Thus, I think we'd only want to apply a filter for abnormally high CV, if at all.  This would remove contigs with the most extreme variation in the mean depth of cvg among positions. 
![image](https://user-images.githubusercontent.com/12803659/66784622-9e541780-eea0-11e9-9725-3fea1528f19b.png)

Digging deeper, the contigs that are isolated by this potential filter tend to have low cvg at most positions, but a minority of positions with higher coverage.  The majority are caught by the min depth of cvg filter.  Those that are not caught by the low cvg fltr have particularly high depths of coverage for most positions, but very low cvg for a few positions at the end of the sequence.  I don't think we want to filter these because they appear to be good contigs with reads that appear to be around 10bp shorter than the average.  e.g.: 

```bash
dDocent_Contig_3689     76      137.143 12245.8
dDocent_Contig_3689     77      136.095 11948.1
dDocent_Contig_3689     78      136.333 11963.8
dDocent_Contig_3689     79      135.905 11991.3
dDocent_Contig_3689     80      136.81  12102.1
dDocent_Contig_3689     81      135.19  11761
dDocent_Contig_3689     82      137.81  12335.4
dDocent_Contig_3689     83      137.143 12132.8
dDocent_Contig_3689     84      137     12105.4
dDocent_Contig_3689     85      136.905 12099.2
dDocent_Contig_3689     86      135.143 11809.6
dDocent_Contig_3689     87      137     12200.9
dDocent_Contig_3689     88      134.333 11601.7
dDocent_Contig_3689     89      136.143 12020.6
dDocent_Contig_3689     90      129.19  10873.7
dDocent_Contig_3689     91      134.714 11635
dDocent_Contig_3689     92      134.429 11691.8
dDocent_Contig_3689     93      135.571 11872.4
dDocent_Contig_3689     94      135.286 11889.3
dDocent_Contig_3689     95      135     11675.6
dDocent_Contig_3689     96      135.81  11905.1
dDocent_Contig_3689     97      134.095 11797.8
dDocent_Contig_3689     98      136.286 12016.6
dDocent_Contig_3689     99      135.333 11915.4
dDocent_Contig_3689     100     133.762 11586.3
dDocent_Contig_3689     101     133     11382.5
dDocent_Contig_3689     102     134.714 11603.3
dDocent_Contig_3689     103     134.857 11783.8
dDocent_Contig_3689     104     134.238 11612.5
dDocent_Contig_3689     105     132.667 11244.8
dDocent_Contig_3689     106     134.81  11868
dDocent_Contig_3689     107     132.619 11414.3
dDocent_Contig_3689     108     131.714 11277.2
dDocent_Contig_3689     109     129.429 10998.5
dDocent_Contig_3689     110     132.952 11565.5
dDocent_Contig_3689     111     132.571 11388.4
dDocent_Contig_3689     112     133.571 11534.7
dDocent_Contig_3689     113     134.476 11753.7
dDocent_Contig_3689     114     132.905 11361.7
dDocent_Contig_3689     115     129.524 10768
dDocent_Contig_3689     116     1.55556 0.527778
dDocent_Contig_3689     117     1.7     0.455556
dDocent_Contig_3689     118     1.63636 0.654545
dDocent_Contig_3689     119     1.92857 1.76374
dDocent_Contig_3689     120     1.85714 1.51648
dDocent_Contig_3689     121     2       1.53846
dDocent_Contig_3689     122     2       2
dDocent_Contig_3689     123     2.14286 2.13187
dDocent_Contig_3689     124     1.91667 0.992424
dDocent_Contig_3689     125     2       1
dDocent_Contig_3689     126     1.53846 0.769231
dDocent_Contig_3689     127     1.55556 0.777778
dDocent_Contig_3689     128     1.75    0.931818
```

## Filtering for extremes in the mean of the CV for each position

This filter would identify  contigs that have extremely high or low CV in the mean depth of coverage by position across individuals.  Here we calculate the mean of the CV in mean cvg by position.  The condition that would make this value particularly high is extreme differences in cvg by individual.  

The mean of the CV in cvg is approximately normally distributed.  Again, I see no reason for a minimum cutoff.  
![image](https://user-images.githubusercontent.com/12803659/66784751-f3902900-eea0-11e9-9f7d-71c204c3af6e.png)

In visually inspecting the most extreme contigs, there tends to be low depth of cvg in most individuals with extremely high depth of a coverage in 1 or a few.  I think it is a good idea to employ this filter with an upper cutoff.  

##  Filtering for extremes in the CV of the CV for each position

Most of the extremes for this filter fail the mean of mean cvg. Those that don't fail the mean cvg filter look similar to the extremes for the cv of the mean of mean cvg.  Thus, I don't think we should filter for this statistic.
![image](https://user-images.githubusercontent.com/12803659/66788461-47a10a80-eead-11e9-9c52-6bdb38a133bd.png)

I pushed a new version of `fltrVCF` to GitHub that allows us to filter by CV of mean coverage (added to fltr041) 
Files that changed: `fltrVCF.bash`, `config.4.all`, `scripts/plotFltr041.R`

Be sure to copy the additional 041 lines into your configs:
```bash
	041 custom bash				0		#Remove contigs with lower CV of mean depth across sites, percentile [0]
	041 custom bash				0.99	#Remove contigs with higher CV of mean depth across sites, percentile [0.99]

```





