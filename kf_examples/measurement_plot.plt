reset
set terminal pngcairo size 800,400 font ",15"

#set cbrange [ 0.0 : 1.0 ] noreverse nowriteback
#set size ratio -1

filename(n) = sprintf("mprof%08d.out", n)
filenamef(n) = sprintf("mprof%08d_filtered.out", n)
fluxfile = "mflux.out"
name(n)=sprintf("img%05d.png",n)

set output "measurements.png"

dt = 1.8e-2

set yrange [0:1]
set y2range [-2.e-5:2.e-5]
set xrange [0:5]
set xlabel "t (h)"

set ytics nomirror tc lt 1
set ylabel "Surface moisture d.b. (/)" tc lt 1

set y2tics nomirror tc lt 2
set y2label "Total flux (kg/m^2s)" tc lt 2

set nokey

j = 1
plot fluxfile u ($2/3600.):6 w l lt 2 lw 2 axis x1y2,\
 "" u ($2/3600.):3 w l lt 1 lw 2 axis x1y1



set output