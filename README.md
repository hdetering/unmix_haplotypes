# Unmix Haplotypes
From a given set of observed haplotypes, infer a minimal set of atomic haplotypes that can serve as their base.

## Getting started
```
git clone https://github.com/hdetering/unmix_haplotypes.git
make
unmix_haplotypes <haplotypes.csv>
```

## Input

Mutation haplotypes are expected to be presented as 0/1 rows in a comma-separated (CSV) file containing row and column IDs.

example:
```
,loc1,loc2,loc3,loc4,loc5
sample_A,1,1,1,0,1
sample_B,0,1,0,1,1
sample_C,0,1,0,1,0
```
