set terminal qt size 1000,550
set xlabel 't' font ',22' offset 0,-1
set ylabel 'S' font ',22' offset -6,0 rotate by 0
set title 'Action' font ',22'
set xtics font ',18'
set ytics font ',18'
set key font ',16'
set key outside
set lmargin 18
set bmargin 5
set grid x,y
set xrange [-1:15]

nn=11
mumax=-3.0

array col1[nn]
array col2[nn]
array m02[nn]
colors="red orange yellow green cyan blue purple brown magenta gray black"

do for [i=1:nn] {
set style line i lc rgb word(colors, i) lw 1.5 pt 1
col1[i]=1+i
col2[i]=nn+1+i
m02[i]=mumax*(i-1.)/(nn-1)
}

plot for [i=1:nn] 'tvssphi.dat' u 1:col1[i]:col2[i] w errorbars title sprintf('μ^2=%.1f',m02[i]) linestyle i
pause -1
