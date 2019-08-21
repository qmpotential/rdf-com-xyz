#!/usr/bin/gnuplot

clear
set terminal pngcairo dashed enhanced size 1024, 768 font 'Times-Roman,18'
set output "SLF_LiBF4_038_019-pos-0.png"
set grid
set xlabel 'r(COM-COM), Å'          font 'Times-Roman,22'
set ylabel 'g(r)'                   font 'Times-Roman,22'
set y2label 'n(r)'                  font 'Times-Roman,22'
set xrange   [  2.00:  8.00] 
set yrange   [  0.00:  4.00] 
set y2range  [  0.00:  4.00] 
set format y "%3.1f"
set format y2 "%3.1f"
set xtics  
set ytics  nomirror  
set y2tics  nomirror  
set mxtics   2.0
set mytics   2.0
set key top left
set pointsize 2
plot \
     \
'SLF_LiBF4_038_019-pos-0.RDF'  u 1:2 w lp lc rgb "blue"   lw 1 lt 1 pt 19 title 'SLF-SLF', \
'SLF_LiBF4_038_019-pos-0.RDF'  u 1:4 w lp lc rgb "black" lw 1 lt 2 pt 19 title 'BF4-BF4', \
'SLF_LiBF4_038_019-pos-0.RDF'  u 1:6 w lp lc rgb "red"  lw 1 lt 3 pt 19 title 'Li-Li', \
'SLF_LiBF4_038_019-pos-0.RDF'  u 1:3 w l lc rgb "blue"   lw 2 lt 1 notitle 'SLF-SLF', \
'SLF_LiBF4_038_019-pos-0.RDF'  u 1:5 w l lc rgb "black" lw 2 lt 2 notitle 'BF4-BF4', \
'SLF_LiBF4_038_019-pos-0.RDF'  u 1:7 w l lc rgb "red"  lw 2 lt 3 notitle 'Li-Li'



