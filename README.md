# HPC_BSG
Performance comparison of an implementation of a variant of the binary search algorithm with gap-encoding and search executions parallelization

For ./binarySearch/binarySearch.cpp
execute as: ./bsearch <n> <prefixResult> <BLK> <NORMAL flag> [<sigma>]
Example:
./bSearch 50 ./ 8 1 10
  
For ./binarySearchParallel/binarySearchP.cpp
execute as: ./bsearch <n> <prefixResult> <BLK> <NORMAL flag> [<sigma>] <threads>
Example:
./bSearch 50 ./ 8 1 10 4
