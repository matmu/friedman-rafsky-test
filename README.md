friedman-rafsky-test
====================
Distribution-based Similarity Measure for Multi-dimensional Point Sets



Introduction
-----------

This is an easy to use Perl implementation of the Friedman Rafsky Test as described in 

_EL Saddik, Abdulmotaleb; Vuong, Son; Griwodz, Carsten; Del Bimbo, Alberto; Candan, K. Selcuk; Jaimes, Alejandro et al. (**2008**): **Distribution-based similarity measures for multi-dimensional point set retrieval applications.** 
In: Proceeding of the 16th ACM international conference on Multimedia - MM '08: ACM Press, S. 429._

The code is written into a single script file called `fr_test.pl` and consists of a workflow with 5 steps:

  * Compute distance matrix,
  * Create minimal spanning tree MST using Prim Algorithm
  * Compute number of runs
  * Compute mean, variance, permutation parameter and quantity W
  * Compute similarity measure (p-Value)

Demo
-----------
Two example datasets _dataset_a.txt_ and _dataset_b.txt_ are located in the same directory as the actual script. An examplary similarity score can be computed via the command
```
./fr_test.pl dataset_a.txt dataset_b.txt 1 -1
```

Perl dependencies
-----------
```
Statistics::Distributions
List::Util qw(sum)
```
Input
-----------
The input takes of two tab-delimited flat files (matrices l x n and m x n) with n attributes, followed by two flags *{-1,1}*, one concerning file headers and the other for printing the distance table.

**Example:**
```
./fr_test.pl dataset1.txt dataset1.txt 1 1
```

Output
-----------
  Similarity score (same as p-Value) is within *[0,1]*. A score of 1 denotes the highest similarity.

**Example:**
```
Total weight of MST: 25528.1492793117
Number of FR-Runs: 38
FR-Permutation parameter C: 130
FR-Variance: 24.5808601478705
FR-Mean: 51
FR-Quantity W: -2.62207321631284
Datasets 'A' and 'B' have a similarity score of 0.436999999999999
```
