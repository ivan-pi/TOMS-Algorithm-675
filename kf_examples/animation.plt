set terminal pngcairo size 800,400 font ",15"
#set size ratio -1

#unset border
#unset xtics
#unset ytics

#set cbrange [ 0.0 : 1.0 ] noreverse nowriteback
#set size ratio -1

filename(n) = sprintf("mprof%08d.out", n)
filenamef(n) = sprintf("mprof%08d_filtered.out", n)
fluxfile = "mflux.out"
name(n)=sprintf("img%05d.png",n)

i=0
istep=100
k=100001
dt = 1.8e-2

load "loop.plt"
