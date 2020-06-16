set terminal postscript eps size 14cm, 14cm
set output "huvetaw.eps"

set xlabel "x, m"
set ylabel "y, m"

set multiplot layout 3,2
  splot "res.out" u 1:2:3 title 'h'
  splot "res.out" u 1:2:4 title 'u'
  splot "res.out" u 1:2:5 title 'v'
  splot "res.out" u 1:2:6 title 'eta'
  splot "res.out" u 1:2:7 title 'w'
unset multiplot
