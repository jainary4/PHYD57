set terminal png size 1400,900 font 'Arial,12'
set output 'convergence_plot.png'

set xlabel 'Number of Samples N' font 'Arial,14'
set ylabel 'Relative Error' font 'Arial,14'
set title 'Monte Carlo Convergence Analysis' font 'Arial,16'

# Use log scale for both axes
set logscale xy

# Grid
set grid

# Key (legend)
set key top right font 'Arial,11'

# Add reference line N^(-0.5) for comparison
f(x) = 0.5 * x**(-0.5)

# Plot the data
plot 'convergence.dat' using 1:2 with linespoints lw 2.5 pt 7 ps 1.5 lc rgb '#e74c3c' title 'P(|t_1-t_2|>1/3)', \
     'convergence.dat' using 1:3 with linespoints lw 2.5 pt 7 ps 1.5 lc rgb '#3498db' title 'E[min(t_1,t_2)]', \
     'convergence.dat' using 1:4 with linespoints lw 2.5 pt 7 ps 1.5 lc rgb '#2ecc71' title 'E[|t_1-t_2|]', \
     'convergence.dat' using 1:5 with linespoints lw 2.5 pt 7 ps 1.5 lc rgb '#f39c12' title 'E[1-max(t_1,t_2)]', \
     f(x) with lines lw 2 lc rgb 'black' dt 3 title 'Reference: N^{-0.5}', \
     1e-5 with lines lw 2 lc rgb 'red' dt 2 title 'Threshold: 10^{-5}'
