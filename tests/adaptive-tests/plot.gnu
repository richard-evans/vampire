set yrange [0:1]
plot 'adaptive/output' u 1:2 t 'adaptive' w lp, \
     'default/output' u 1:2 t 'default' w lp
