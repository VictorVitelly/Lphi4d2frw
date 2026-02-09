set terminal qt size 1000,600
set xlabel 'μ^2' font ',22' offset 0,-1
set ylabel 'S' font ',22' offset -3,0 rotate by 0
set title 'Action' font ',20'
set xtics font ',18'
set ytics font ',18'
set key font ',16'
set grid x,y
#set key outside
set lmargin 15
set bmargin 5
set key outside
set xrange [-3.25:0.25]

nn=10

colors="red orange yellow green cyan blue purple brown magenta gray black"

array col1[nn]
array col2[nn]
array tt[nn]

do for [i=1:nn] {
set style line i lc rgb word(colors, i) lw 1.8 pt 1 
col1[i]=2.*i
col2[i]=1+2.*i
tt[i]=i
}

tt[1]=0
tt[2]=1
tt[3]=3
tt[4]=5
tt[5]=7
tt[6]=9
tt[7]=11
tt[8]=13
tt[9]=14
tt[10]=15

plot for [i=1:9] 'm02vssphi.dat' u 1:col1[i]:col2[i] w errorbars title sprintf('t=%.0f',tt[i]) linestyle i, 'action.dat' w errorbars title 'Full lattice' pt 2 lc rgb "black", 'actionstatic.dat' w lines title 'Static 2d' lc rgb "black", 'actionstatic1d.dat' w lines title 'Static 1d'
pause -1
