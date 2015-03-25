# Simple2.0
Scripts adapted from http://simple.stanford.edu/ to run on HPC with SLURM system

##### Example commands to run in parallel:

```distr_simple.pl prg=simple_rndrec stk=STKrib.spi box=120 smpd=3.54 nthr=16 npart=2 queue=grant```

```distr_simple_prime.pl stk=STKrib.spi box=120 smpd=3.54 nthr=16 npart=2 queue=grant vol1=startvol_state1.spi dynlp=yes oritab=rndoris.txt ring2=50```
