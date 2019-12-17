I'm assuming that you are in your `mkBAM` dir.

1. Edit the `dDocentHPC.sbatch` file as follows
```bash
#bash dDocentHPC.bash mkBAM config.4.all
bash dDocentHPC.bash fltrBAM config.4.all
```
2. run the filtering
```bash
sbatch dDocentHPC.sbatch
```
I'm noticing better performance in mapping and filtered mapping with the contigs that don't have `NNNNNNNNNNNNNNNNNNNN`, i.e., the contigs that assembled.  90% of my contigs assembled and did not have the `NNNNNNNNNNNNN`
