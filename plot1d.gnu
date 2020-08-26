set terminal postscript eps color size 14cm, 14cm
set output "huvetaw1d.eps"
set style data lines
#set pointsize 0.5

set xlabel "x, m"

set multiplot layout 3,2
  set ylabel "h, m"#; set yrange[0.9:1.9]
  plot "res.out" u 2:3 title 'h'

  set ylabel "u, m/s"#; set yrange[-1:1]
  plot "res.out" u 2:4 title 'u'

  set ylabel "v, m/s"
  plot "res.out" u 2:5 title 'v'

  set ylabel "eta, m"#; set yrange[0.9:1.9]
  plot "res.out" u 2:6 title 'eta'

  set ylabel "w, m/s"#; set yrange[-1:1]
  plot "res.out" u 2:7 title 'w'
unset multiplot
