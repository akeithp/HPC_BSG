set terminal postscript eps color solid "Helvetica" 16
#set terminal postscript eps color colortext
set output "./bSearchP.eps"
set title "Array length vs. Query time Parallel"
set xlabel "n"
set ylabel "Time (microsecs)"
set style func linespoints
set pointsize 1
set key left
plot [] [] \
	'./resultsP/bSearchBLK1000' using 1:4 title "BS BLK1000" with linespoints lt 0 pt 7 lw 0.5, \
	'./resultsP/bSearchScanBLK1000' using 1:4 title "BSSCN BLK1000" with linespoints lt 1 pt 8 lw 0.5,\
	'./resultsP/bSearchGapBLK1000' using 1:4 title "BSGSCN BLK1000" with linespoints lt 2 pt 9 lw 0.5,\
	'./resultsP/bSearchBLK20000' using 1:4 title "BS BLK20000" with linespoints lt 3 pt 5 lw 0.5, \
	'./resultsP/bSearchScanBLK20000' using 1:4 title "BSSCN BLK20000" with linespoints lt 4 pt 6 lw 0.5, \
	'./resultsP/bSearchGapBLK20000' using 1:4 title "BSGSCN BLK20000" with linespoints lt 5 pt 4 lw 0.5

