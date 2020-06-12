j = i/istep
set output name(j)
#print filename(i)
#print name(j)
#set title sprintf("t = %4.2f h",i*dt/3600.)
#unset colorbox

set multiplot layout 1,2 columnsfirst
set yrange [0:0.5]
set xrange [0:0.5]

set xlabel "x (mm)"
set ylabel "Moisture"
plot filename(i) using ($1*1000):2 with l lw 2 notitle

set yrange [0:1]
set xrange [0:5]
set xlabel "t (h)"
set ylabel "Moisture"
set arrow 1 from i*dt/3600.,0 to i*dt/3600.,1 ls 1 lc rgb "black" nohead
plot fluxfile every ::::j u ($2/3600.):3 w l lw 2 notitle,\
 "" every ::::j u ($2/3600.):4 w l lw 2 notitle,\
 "" every ::::j u ($2/3600.):5 w l lw 2 notitle,\
 "" u ($2/3600.):7 w l dt 2 notitle,\
 "" every ::::j u ($2/3600.):9 w l lw 2 notitle

unset arrow 1
unset multiplot
replot

i=i+istep
if (i < k) reread
