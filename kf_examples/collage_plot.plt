reset
set terminal pngcairo size 800,400 font ",15"

#set cbrange [ 0.0 : 1.0 ] noreverse nowriteback
#set size ratio -1

filename(n) = sprintf("mprof%08d.out", n)
filenamef(n) = sprintf("mprof%08d_filtered.out", n)
fluxfile = "mflux.out"
name(n)=sprintf("img%05d.png",n)

set output "collage.png"

i=0
i1=50000
i2=400000
i3=850000
i4=950000

dt = 1.8e-2

set multiplot layout 1,2 columnsfirst
set yrange [0:0.5]
set xrange [0:0.5]

set xlabel "x (mm)"
set ylabel "Moisture"
plot filename(i1) using ($1*1000):2 with l lw 2 notitle,\
     filename(i2) using ($1*1000):2 with l lw 2 notitle,\
     filename(i3) using ($1*1000):2 with l lw 2 notitle,\
     filename(i4) using ($1*1000):2 with l lw 2 notitle,\

set yrange [0:1]
set xrange [0:5]
set xlabel "t (h)"

unset ylabel
set arrow 1 from i1*dt/3600.,0 to i1*dt/3600.,1 dt 2 lw 2 lc 1 nohead
set arrow 2 from i2*dt/3600.,0 to i2*dt/3600.,1 dt 2 lw 2 lc 2 nohead
set arrow 3 from i3*dt/3600.,0 to i3*dt/3600.,1 dt 2 lw 2 lc 3 nohead
set arrow 4 from i4*dt/3600.,0 to i4*dt/3600.,1 dt 2 lw 2 lc 4 nohead

set key at 1,1.5 font ",11" box opaque
set border back

j = 1
plot fluxfile u ($2/3600.):7 w l lw 2 lc rgb "black" title "RH (%)",\
 "" u ($2/3600.):9 w l lw 2 lc rgb "blue" title "Average moisture"


unset multiplot

set output
