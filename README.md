# pairwise-openmp
Parallel pairwise alignment using OpenMP using various data partition scheme:
- columnwise
- rowwise
- antidiagonal
- blocked antidiagonal

For performance comparison dan results, *see* the [paper](paper.pdf).

## Compile and Run Test

```
gcc columnwise.c -o columnwise -fopenmp
./columnwise < dna-16.in
```

## Debug Mode

Compile with `-D` option.

```
gcc -DDEBUG columnwise.c -o columnwise -fopenmp
```

Run with any two strings with specified length.

```
./columnwise
6
4
auriza
auzi
```

Output.

```
auriza
au__zi

  0  -3  -6  -9 -12 
 -3   1  -2  -5  -8 
 -6  -2   2  -1  -4 
 -9  -5  -1   1  -2 
-12  -8  -4  -2   2 
-15 -11  -7  -3  -1 
-18 -14 -10  -6  -4 

  a u z i
a \ _ _ _
u | \ _ _
r | | \ \
i | | \ \
z | | \ |
a \ | | \
Time: 0.000402 s	Score: -4
```
