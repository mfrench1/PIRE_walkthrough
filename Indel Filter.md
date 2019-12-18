Hi All,

I don't have much time but I did put together this line of code which will result in a list of contigs that 2 or more indels of 3 or more nucleotides.  Actually, I made it so that it's a list of regex patterns for the aforementioned contigs that can be passed to grep.

Note that I applied this to the results of fltr30 because that zeroed in on the positions we will be targeting with the probes.  You should do the same

First, make sure it works for you because this might not be bullet proof.  I rushed it.

```bash
grep 'TYPE=ins\|TYPE=del' PIRE_SiganusSpinus.L.5.5.Fltr30.1.vcf | cut -f1,8 | tr ";" "\t" | cut -f1,15 | sed  -e 's/,[12],/,/g' -e 's/,[12],/,/g' -e 's/,[12]$//g' -e 's/,[12]$//g' -e 's/=[12],/=/g' -e 's/=[12],/=/g' -e 's/=[12]$/=/g' | grep -v 'LEN=$' | cut -f1 | uniq -c | tr -s " " "\t" | cut -f2-3 | grep -vP '^1\t' | cut -f2 | less

#now search inside less for the following patterns independently:  =1$  =2$ =1, =2,  ,1,  ,2,
# you should not find any lines that match.  remember to goto top of doc after each search using "g"
```

It it looks good, then:

```bash
grep 'TYPE=ins\|TYPE=del' PIRE_SiganusSpinus.L.5.5.Fltr30.1.vcf | cut -f1,8 | tr ";" "\t" | cut -f1,15 | sed  -e 's/,[12],/,/g' -e 's/,[12],/,/g' -e 's/,[12]$//g' -e 's/,[12]$//g' -e 's/=[12],/=/g' -e 's/=[12],/=/g' -e 's/=[12]$/=/g' | grep -v 'LEN=$' | cut -f1 | uniq -c | tr -s " " "\t" | cut -f2-3 | grep -vP '^1\t' | cut -f2 |  sed -e 's/^/\^>/' -e 's/$/\t/' > PIRE_SiganusSpinus.L.5.5.Fltr33.remove.contigs 
```

You can then use grep on the reference genome with this file of search patterns to remove the contigs with too many indels

```bash
cat PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.371-618bp.0-8TGCAGG.fasta | paste - - | grep -vf PIRE_SiganusSpinus.L.5.5.Fltr33.remove.contigs | tr "\t" "\n" > PIRE_SiganusSpinus.L.5.5.probes4development2.noMSATS.noNNNN.371-618bp.0-8TGCAGG.lessthan2indelsof3nt.fasta
```

Can somebody take a crack at modifying the above code to remove contigs with more than X indels ( period )?  Perhaps 3 or more indels?

Perhaps we should also try to remove loci that have an indel larger than Y. I sent an inquiry to Arbor for advice on the thresholds.
