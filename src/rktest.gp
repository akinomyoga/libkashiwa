#!/usr/bin/gnuplot
# -*- mode:sh -*-

set terminal pdf color enhanced size 6,8 #4.8,6.4
set output "out/rk/rktest.pdf"
set title 'rktest: Eq.: {/Times ~{/Times-Italic x}{.9.} = tan {/Times-Italic x} + 1}, IC: {/Times {/Times-Italic x}(0) = 1}.'
#set title 'rktest: Eq.: {/Times ~{/Times-Italic x}{.9.} = 1/{/Times-Italic x}}, IC: {/Times {/Times-Italic x}(0) = 1}.'
set log x
set log y
set format x "10^{%L}"
set format y "10^{%L}"
set xlabel '{/Times {/Times-Italic x}(0.1)} abs(error)'
set ylabel '#{/Times-Italic f} eval'
set xrange [] reverse
set key left top reverse Left
yaxis(y,e) = y
#yaxis(x,e) = x*abs(e)**(1.0/3.0)
#yaxis(x,e) = x*abs(e)**(1.0/6.0); set yrange [0.2:40] # 6次法で規格化した計算コスト
plot \
  "out/rk/rkeuler.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Euler', \
  "out/rk/rkmid.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'midpoint', \
  "out/rk/rkral.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Ralston', \
  "out/rk/rkheun.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Heun', \
  "out/rk/runge3.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Runge3', \
  "out/rk/rkheun3.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Heun3', \
  "out/rk/rkral3.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Ralston3', \
  "out/rk/rkk3.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Kutta3', \
  "out/rk/rkrk4.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'RK4', \
  "out/rk/rkgill.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Gill', \
  "out/rk/rkk38.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Kutta 3/8', \
  "out/rk/b5v3.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v3', \
  "out/rk/b5v2.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v2', \
  "out/rk/b5v1.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v1', \
  "out/rk/hammud6.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Hammud6', \
  "out/rk/shanks7.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Shanks7', \
  "out/rk/rkcv7.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Cooper Verner 7', \
  "out/rk/dop853.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'DOP853 fixed h', \
  "out/rk/mdop853.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'DOP853', \
  "out/rk/rkcv8.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Cooper Verner 8'

yaxis(fe,err) = (-log(abs(err)/0.1)/log(abs(fe))); set nolog y; unset format y; set ylabel 'Score'; set yrange [0:*] reverse; set xrange [*:10]; set key left bottom
plot \
  "out/rk/rkeuler.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Euler', \
  "out/rk/rkmid.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'midpoint', \
  "out/rk/rkral.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Ralston', \
  "out/rk/rkheun.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Heun', \
  "out/rk/runge3.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Runge3', \
  "out/rk/rkheun3.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Heun3', \
  "out/rk/rkral3.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Ralston3', \
  "out/rk/rkk3.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Kutta3', \
  "out/rk/rkrk4.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'RK4', \
  "out/rk/rkgill.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Gill', \
  "out/rk/rkk38.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Kutta 3/8', \
  "out/rk/b5v3.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v3', \
  "out/rk/b5v2.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v2', \
  "out/rk/b5v1.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v1', \
  "out/rk/hammud6.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Hammud6', \
  "out/rk/shanks7.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Shanks7', \
  "out/rk/rkcv7.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Cooper Verner 7', \
  "out/rk/dop853.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'DOP853 fixed h', \
  "out/rk/mdop853.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'DOP853', \
  "out/rk/rkcv8.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Cooper Verner 8'
