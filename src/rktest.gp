#!/usr/bin/gnuplot
# -*- mode:sh -*-

outdir="../out/rk"

set terminal pdf color enhanced size 6,8 #4.8,6.4
set output outdir."/rktest.pdf"
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
  outdir."/rkeuler.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Euler', \
  outdir."/rkmid.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'midpoint', \
  outdir."/rkral.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Ralston', \
  outdir."/rkheun.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Heun', \
  outdir."/runge3.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Runge3', \
  outdir."/rkheun3.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Heun3', \
  outdir."/rkral3.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Ralston3', \
  outdir."/rkk3.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Kutta3', \
  outdir."/rkrk4.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'RK4', \
  outdir."/rkgill.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Gill', \
  outdir."/rkk38.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Kutta 3/8', \
  outdir."/b5v3.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v3', \
  outdir."/b5v2.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v2', \
  outdir."/b5v1.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher5v1', \
  outdir."/hammud6.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Hammud6', \
  outdir."/butcher6.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Butcher6', \
  outdir."/shanks7.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Shanks7', \
  outdir."/rkcv7.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Cooper Verner 7', \
  outdir."/dop853.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title '(DOP853)', \
  outdir."/rkcv8.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title 'Cooper Verner 8', \
  outdir."/verner9.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 title '(Verner9)'

yaxis(fe,err) = (-log(abs(err)/0.1)/log(abs(fe)))
set nolog y; unset format y; set ylabel 'Score'; set yrange [0:*] reverse; set xrange [*:10]; set key left bottom
plot \
  outdir."/rkeuler.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#222222' title 'Euler', \
  outdir."/rkmid.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#884400' title 'midpoint', \
  outdir."/rkral.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#CCAAAA' title 'Ralston', \
  outdir."/rkheun.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#880000' title 'Heun', \
  outdir."/runge3.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#AACCAA' title 'Runge3', \
  outdir."/rkheun3.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#008844' title 'Heun3', \
  outdir."/rkral3.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#AACCCC' title 'Ralston3', \
  outdir."/rkk3.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#006600' title 'Kutta3', \
  outdir."/tvdrk3.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#338800' title 'TVD-RK3', \
  outdir."/tvdrk43.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#338800' dt (16,4) title 'TVD-RK3(4)', \
  outdir."/rkrk4.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#000088' title 'RK4', \
  outdir."/rkgill.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#6666AA' title 'Gill', \
  outdir."/rkk38.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#0000FF' title 'Kutta 3/8', \
  outdir."/b5v3.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#880088' title 'Butcher5v3', \
  outdir."/b5v2.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#CCAACC' title 'Butcher5v2', \
  outdir."/b5v1.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#FF00FF' title 'Butcher5v1'

yaxis(fe,err) = (-log(abs(err)/10.0)/log(abs(fe)))
set ylabel 'Score2'
set yrange [3:*] reverse
plot \
  outdir."/b5v1.txt"     u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#FF00FF' title 'Butcher5v1', \
  outdir."/hammud6.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#CCCC00' title 'Hammud6', \
  outdir."/butcher6.txt" u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#CC8800' title 'Butcher6', \
  outdir."/shanks7.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#88AACC' title 'Shanks7', \
  outdir."/rkcv7.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#00BBBB' title 'Cooper Verner 7', \
  outdir."/dop853.txt"   u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#666666' title '(DOP853)', \
  outdir."/mdop853.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#666666' dt (4,4) title 'DOP853', \
  outdir."/rkcv8.txt"    u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#000000' title 'Cooper Verner 8', \
  outdir."/verner9.txt"  u (abs($5)):(yaxis($1,$5)) w lp ps 0.5 lc rgb '#FF0000' title '(Verner9)'
