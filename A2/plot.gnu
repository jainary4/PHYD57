set terminal png size 1400,900 font 'Arial,12'
set output 'radial_velocity.png'
set xlabel 'Mean Anomaly nt (radians)' font 'Arial,14'
set ylabel 'Radial Velocity V_r (km/s)' font 'Arial,14'
set title 'Stellar Radial Velocity Due to Exoplanet (M=1193399999999999912998055444480.0 M_{sun}, a=1.32 AU, e=0.68)' font 'Arial,16'
set grid
set key top right font 'Arial,11'
set xrange [0:10]
plot 'omega0.dat' using 1:2 with lines lw 2.5 lc rgb '#e74c3c' title 'omega = 0 rad', \
     'omega1.dat' using 1:2 with lines lw 2.5 lc rgb '#3498db' title 'omega = 1 rad', \
     'omega2.dat' using 1:2 with lines lw 2.5 lc rgb '#2ecc71' title 'omega = 2 rad'
