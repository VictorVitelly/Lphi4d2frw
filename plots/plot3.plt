set terminal qt size 880,550
set xlabel 'μ^2' font ',22' offset 0,-1
set ylabel 'χ_t' font ',22' offset -3,0 rotate by 0
set title 'Robin Boundary Conditions' font ',20'
set xtics font ',18'
set ytics font ',18'
set key font ',18'
set grid x,y
#set key outside
set lmargin 15
set bmargin 5
set xrange [-3.25:0.25]

nn=10

colors="red orange yellow green cyan blue purple brown magenta gray black"

array col1[nn]
array col2[nn]
array tt[nn]

do for [i=1:nn] {
set style line i lc rgb word(colors, i) lw 1.8 pt 1 
col1[i]=2.*i
col2[i]=1+2*i
tt[i]=i
}

plot for [i=1:nn] 'susceptibility.dat' u 1:col1[i]:col2[i] w errorbars title sprintf('t=%.0f',tt[i]) linestyle i 
pause -1
