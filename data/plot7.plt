set title 'Correlation length' font ',20'
#set title 'Renormalized mass' font ',20'
set xlabel 'μ^2' font ',20' offset 0,-1
set ylabel 'ξ' font ',20' offset -3,0 rotate by 0
set xtics font ',18'
set ytics font ',18'
set key font ',16'
set grid x,y
#set key outside
set lmargin 15
set bmargin 5
#set key left
set xrange [-3.25:0.25]

p 'longt0.dat' u 1:2:3 w errorbars title 't=0', 'longt7.dat' u 1:2:3 w errorbars title 't=7', 'longt15.dat' u 1:2:3 w errorbars title 't=15'

#p 'mrt0.dat' u 1:2:3 w errorbars title 't=0', 'mrt7.dat' u 1:2:3 w errorbars title 't=7', 'mrt15.dat' u 1:2:3 w errorbars title 't=15'

pause -1
