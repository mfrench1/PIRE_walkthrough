I have noticed a few abnormalities in the reference genome and we should probably address them prior to probe development by filtering some of the contigs.  This can be done either before or after the VCF filter is applied to the reference genome.

## Microsatellites

Approximately 1/30th of my ref genome has microsatellite motifs.  I'm concerned that the indels associated with microsats would promote allelic dropout.  I've written about how to apply this filter at the end of filtering the VCF file.  

## Contigs with NNNNNNNNN*

The vast majority of my contigs were assembled completely, have no NNNNNNNN and they're around 500-600bp.  However, my contigs with NNNNNNNN are much larger, on average (>800bp).  Basically, there is no reason to expect this behavior and I think that something is going wrong in the way dDocent interfaces with rainbow on these types of assemblies.  Jon and I identified a problem with this part of the code earlier in 2019 and it has changed.  He has additionally changed the way single digest data is handled and I suspect that he doesn't work with it as much as double digest.

200/30000 contigs have NNNNNNN and I think we should filter them from consideration. If others have a high proportion of contigs with NNNNNNN we can revisit.

```bash
# convert fasta to tab delim
cat PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.fasta | paste - - > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.tsv

#filter contigs with NNNNNN and convert back to fasta
grep -v 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.tsv | tr "\t" "\n" > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta
```

## Contig Length

Illumina sequencers need short sequences, and if we are assembling a contig with >800 bp, it's potentially suspect.  I think we should apply a contig length filter.  Either a percentile based one or a hard cutoff around 800 bp. After removing the contigs with NNNNNN, there are very few contigs anywhere near 800 bp.


```bash
# make tsv with name and num bp of each contig
paste <(grep '^>'  PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta) \
	<(grep -v '^>' PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta | awk '{ print length }') | \
	sort -nk2 > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv

# experiment with the contig lengths returned by different percentiles
cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=0.99 '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}'
618
cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=0.999 '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}'
745
cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=0.001 '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}'
328

# apply filter
THRESHOLD=0.01
THRESHOLDb=0.99

THRESHOLD=$(cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=$THRESHOLD '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
THRESHOLDb=$(cut -f2 PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | awk -v PCT=$THRESHOLDb '{all[NR] = $0 } END{print all[int(NR*PCT - 0.5)]}')
awk -v LENb=$THRESHOLDb -v LEN=$THRESHOLD '$2 > LENb || $2 < LEN {print $1;}' PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.tsv | sed -e 's/^/\^/' -e 's/$/\t/' > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.remove.contigs
cat PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.fasta | paste - - > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.tsv
grep -vf PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.lengths.remove.contigs PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.tsv | tr "\t" "\n" > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.${THRESHOLD}-${THRESHOLDb}bp.fasta


```



