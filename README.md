##TrioMDR:
detecting SNP interactions in trio families with model-based multifactor dimensionality reduction

Description :This page includes the codes of TrioMDR which are implemented and tested on R version 3.2.4 by Mr. Liu (Jie Liu, Email:jiel@email.swu.edu.cn) and free for academic usage. You can run it at your own disk. For other purposes, please contact the corresponding Author Dr. Wang (JunWang, Email: kingjun@swu.edu.cn). If you have any problem on using these codes,just contact Dr. Wang via Email.	


Usage:
source("path/TrioMDR.R")
results <- TrioMDR(Trios.file, K=2, classm='PDT', adj.main='FALSE', nperm=5)


Parameters:
Trios.file: The PLINK format ped file. Not that the ped file should only contain numbers. See test.ped as an example (loucs 2 and 7 are disease loci).
classm: classification rule to define H/L, default geno-PDT statistic
adj.main: logical value, adjust marginal effect or not, default FALSE.
nperm: number of permutation time to estimate non-central parameters, default 5


Column Number Description of Trios.file:
1 Pedigree number (must be an integer)
2 Individual ID number (must be an integer)
3 ID of father (0 if this individual is a founder)
4 ID of mother (0 if this individual is a founder)
5 Sex (1=male, 2=female, 0=unknown)
6 Disease status (2=affected, 1=unaffected, 0=unknown)
7 The 1rst allele pair (must be pairs of integers)
8
... ...
n-1 The 1rst allele of the last pair
n


