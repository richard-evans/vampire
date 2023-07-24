set terminal pngcairo
set xlabel 'Input T(K)'
set ylabel 'Spin temperature'

set output 'plot_spin_temp.png'
f(x)=x
p 'output' u 1:4 w lp lw 2 ps 1.5 t 'Instantaneous Ts', '' u 1:5 w lp lw 2 ps 1.5 pt 7  lc rgb 'red' t 'Mean Ts', f(x) w l lw 1.5 lc -1