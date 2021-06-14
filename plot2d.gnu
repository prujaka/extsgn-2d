set terminal postscript eps enhanced color size 14cm, 14cm
set output "huvetaw2d.eps"
set pointsize 0.2

set xlabel "x, m"
set ylabel "y, m"

set multiplot layout 3,2
  splot "res.dat" u 1:2:3 title 'h'
  splot "res.dat" u 1:2:4 title 'u'
  splot "res.dat" u 1:2:5 title 'v'
  splot "res.dat" u 1:2:6 title 'eta'
  splot "res.dat" u 1:2:7 title 'w'
unset multiplot
