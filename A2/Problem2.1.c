#include <stdio.h>
#include<math.h>
#include <stdlib.h>
#define M_SUN (1.989e30)
#define AU_M (1.496e11)
#define SQRT_GM_SUN_PER_AU 29.78 
// Planet data parameters
double M_star;
double M_planet;
double M_total=0.6*M_SUN;
double G = 6.67430e-11;
double e=0.68;
double a= 1.32;
int MAX_ITER=50;
double m_p_over_M = 0.001;
double GM_SUN_AU=29.78;



double solve_kepler(double n, double t) {
    // solves the kepler's equation to find the eccentric anomaly
    double E = n*t;  // Initial guess
    double E_prev;
    int iter = 0;
    
    while(iter<=MAX_ITER){
        E_prev=E;
        E= n*t+e*sin(E_prev);
        iter++;
    }
    return E;
}

double eccentric_to_true_anomaly(double E) {
    // converts the eccentric anomaly to the true anomaly
    double tan_f_2 = sqrt((1.0 + e) / (1.0 - e)) * tan(E / 2.0);
    return 2.0 * atan(tan_f_2);
}

double calc_radius_1(double f, double a) {
    // calculates the radius w.r.t to f
    return a * (1.0 - e * e) / (1.0 + e * cos(f));
}

double calc_radius_2(double E, double a) {
    return a * (1.0 - e * cos(E));
}

double calc_radial_velocity(double t,double n){
    /*This function calculates the radial velocity of a planet orbiting in a elliptical orbit around its star
    the equation is given by e*sin(f)(sqrt(G*M_total/a(1-e^2)))
    we need to first find f to calculate the radial velocity*/
    double E=solve_kepler(n,t);
    double f=eccentric_to_true_anomaly(E);
    double radial_velocity=e*sin(f)*sqrt(G*M_total/(a*(1.0-e*e)));
    printf("The radial velocity of the planet is %f km/s\n", radial_velocity);

    return m_p_over_M*radial_velocity;
}

double calc_perpendicular_velocity(double t,double n){
    /*This function calculates the perpendicular velocity of a planet orbiting in a elliptical orbit around its star
    the equation is given by (1+e*cos(f))*sqrt(G*M_total/a(1-e^2))
    we need to first find f to calculate the perpendicular velocity*/
    double E=solve_kepler(n,t);
    double f=eccentric_to_true_anomaly(E);
    double perpendicular_velocity=(1.0+e*cos(f))*sqrt(G*M_total/(a*(1.0-e*e)));
    printf("The perpendicular velocity of the planet is %f km/s\n", perpendicular_velocity);

    return m_p_over_M*perpendicular_velocity;

}

double radial_velocity_of_star_relative_to_earth(double t, double w, double n){
    /*This function calculates the radial velocity of the star relative to the earth*/
    double E=solve_kepler(n,t);
    double f=eccentric_to_true_anomaly(E);
    double velocity_scale = SQRT_GM_SUN_PER_AU * sqrt(0.6);
    double radial_velocity_of_star_relative_to_earth= velocity_scale*sqrt(1/(a*(1.0-e*e)))*(sin(w+f)*(1+e*cos(f))-e*sin(f)*cos(w+f));
    printf("The radial velocity of the star relative to the earth is %f km/s\n", radial_velocity_of_star_relative_to_earth);

    return -1*(m_p_over_M*radial_velocity_of_star_relative_to_earth);
}

int main(){
double omega_values[3] = {0.0, 1.0, 2.0};
const char* filenames[3] = {"omega0.dat", "omega1.dat", "omega2.dat"};

// Store amplitudes for each omega
double amplitudes[3];
double n = 1.0;
int N=500;

for (int w = 0; w < 3; w++) {
    // opens the file for writing
    FILE *fp = fopen(filenames[w], "w");

    if (fp == NULL) {
        fprintf(stderr, "Error opening file %s\n", filenames[w]);
        return 1;
    }
    
    fprintf(fp, "# Mean Anomaly (rad)    Radial Velocity (km/s)\n");
    fprintf(fp, "# omega = %.1f rad\n", omega_values[w]);
    
    // Calculate for nt = 0 to 10 radians
    double min_vel = 1e10, max_vel = -1e10;
    
    for (int i = 0; i < N; i++) {
        double t = (i / (double)(N - 1)) * 10.0;
        double Vr = radial_velocity_of_star_relative_to_earth(t, omega_values[w], n);
        fprintf(fp, "%.6f    %.6f\n", n * t, Vr);
        
        // Track min and max
        if (Vr < min_vel) min_vel = Vr;
        if (Vr > max_vel) max_vel = Vr;
    }
    
    amplitudes[w] = max_vel - min_vel;
    
    fclose(fp);

}

printf("Generated files: omega0.dat, omega1.dat, omega2.dat\n");
    
    // Create gnuplot script
    FILE *gp = fopen("plot.gnu", "w");
    if (gp != NULL) {
        fprintf(gp, "set terminal png size 1400,900 font 'Arial,12'\n");
        fprintf(gp, "set output 'radial_velocity.png'\n");
        fprintf(gp, "set xlabel 'Mean Anomaly nt (radians)' font 'Arial,14'\n");
        fprintf(gp, "set ylabel 'Radial Velocity V_r (km/s)' font 'Arial,14'\n");
        fprintf(gp, "set title 'Stellar Radial Velocity Due to Exoplanet (M=%.1f M_{sun}, a=%.2f AU, e=%.2f)' font 'Arial,16'\n", 
                M_total, a , e);
        fprintf(gp, "set grid\n");
        fprintf(gp, "set key top right font 'Arial,11'\n");
        fprintf(gp, "set xrange [0:10]\n");
        fprintf(gp, "plot 'omega0.dat' using 1:2 with lines lw 2.5 lc rgb '#e74c3c' title 'omega = 0 rad', \\\n");
        fprintf(gp, "     'omega1.dat' using 1:2 with lines lw 2.5 lc rgb '#3498db' title 'omega = 1 rad', \\\n");
        fprintf(gp, "     'omega2.dat' using 1:2 with lines lw 2.5 lc rgb '#2ecc71' title 'omega = 2 rad'\n");
        fclose(gp);
        printf("Generated: plot.gnu\n");
        printf("Run 'gnuplot plot.gnu' to create radial_velocity.png\n");
    }
    

    
    return 0;
}