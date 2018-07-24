max=8
for i in {1..8}
do
	./binarySearchParallel/bSearch 100000000 ./resultsPN/ 1000 1 10000 $i
	./binarySearchParallel/bSearch 100000000 ./resultsPN/ 20000 1 10000 $i
	./binarySearchParallel/bSearch 100000000 ./resultsPU/ 1000 0 10000 $i
	./binarySearchParallel/bSearch 100000000 ./resultsPU/ 20000 0 10000 $i
done



