#include<stdio.h>
#include<math.h>

// Coulomb's law constant
double k = 9e9;
// Charge of an electron/positron
double e = 1.60217663e-19;
// Distance between the charges (100 fm)
double d = 100e-15;
// Masses
double m_e = 9.1093837e-31;   
double m_p = 1.6726219e-27;   

// Structure to hold particle state
typedef struct {
    double x, y;      // position
    double vx, vy;    // velocity
    double mass;      // mass
    double charge;    // charge
} Particle;

// Function to calculate acceleration on particle i due to all other particles
void calculate_acceleration(Particle particles[], int i, int n_particles, double *ax, double *ay) {
    *ax = 0.0;
    *ay = 0.0;
    
    for (int j = 0; j < n_particles; j++) {
        if (i == j) continue;  
        
        // Distance between particles i and j
        double dx = particles[j].x - particles[i].x;
        double dy = particles[j].y - particles[i].y;
        double r = sqrt(dx*dx + dy*dy);
        
        
        double F_mag = k * particles[i].charge * particles[j].charge / (r * r);
        
        // Force components (repulsive, so force on i points away from j)
        double Fx = -F_mag * dx / r;  
        double Fy = -F_mag * dy / r;
        
        // Acceleration components: a = F/m
        *ax += Fx / particles[i].mass;
        *ay += Fy / particles[i].mass;
    }
}

// RK2 integration step for all 16 variables (4 particles × 4 variables each)
void rk2_step(Particle particles[], int n_particles, double dt) {
    // Arrays to store k1 and k2 values for each particle
    double k1_x[4], k1_y[4], k1_vx[4], k1_vy[4];
    double k2_x[4], k2_y[4], k2_vx[4], k2_vy[4];
    Particle temp_particles[4];
    
    
    for (int i = 0; i < n_particles; i++) {
        double ax, ay;
        calculate_acceleration(particles, i, n_particles, &ax, &ay);
        
        k1_x[i] = particles[i].vx;     // dx/dt = vx
        k1_y[i] = particles[i].vy;     // dy/dt = vy
        k1_vx[i] = ax;                  // dvx/dt = ax = Fx/m
        k1_vy[i] = ay;                  // dvy/dt = ay = Fy/m
    }
    
    
    for (int i = 0; i < n_particles; i++) {
        temp_particles[i] = particles[i]; 
        temp_particles[i].x = particles[i].x + k1_x[i] * dt / 2.0;
        temp_particles[i].y = particles[i].y + k1_y[i] * dt / 2.0;
        temp_particles[i].vx = particles[i].vx + k1_vx[i] * dt / 2.0;
        temp_particles[i].vy = particles[i].vy + k1_vy[i] * dt / 2.0;
    }
    
    
    for (int i = 0; i < n_particles; i++) {
        double ax, ay;
        calculate_acceleration(temp_particles, i, n_particles, &ax, &ay);
        
        k2_x[i] = temp_particles[i].vx;
        k2_y[i] = temp_particles[i].vy;
        k2_vx[i] = ax;
        k2_vy[i] = ay;
    }
    
    
    for (int i = 0; i < n_particles; i++) {
        particles[i].x += k2_x[i] * dt;
        particles[i].y += k2_y[i] * dt;
        particles[i].vx += k2_vx[i] * dt;
        particles[i].vy += k2_vy[i] * dt;
    }
}

// Calculate total kinetic energy
double calculate_kinetic_energy(Particle particles[], int n_particles) {
    double KE = 0.0;
    for (int i = 0; i < n_particles; i++) {
        double v_squared = particles[i].vx * particles[i].vx + 
                          particles[i].vy * particles[i].vy;
        KE += 0.5 * particles[i].mass * v_squared;
    }
    return KE;
}

// Calculate total potential energy 
double calculate_potential_energy(Particle particles[], int n_particles) {
    double PE = 0.0;
    
    // Sum over all unique pairs (i < j)
    for (int i = 0; i < n_particles; i++) {
        for (int j = i + 1; j < n_particles; j++) {
            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;
            double r = sqrt(dx*dx + dy*dy);
            
            
            PE += k * particles[i].charge * particles[j].charge / r;
        }
    }
    return PE;
}

// Calculate initial potential energy for verification
double calc_initial_potential_energy() {
    // E0 = k*e²/d * (4 + sqrt(2))
    // 4 adjacent pairs at distance d
    // 2 diagonal pairs at distance d*sqrt(2)
    double U_adjacent = 4.0 * k * e * e / d;
    double U_diagonal = 2.0 * k * e * e / (d * sqrt(2.0));
    double U_total = U_adjacent + U_diagonal;
    
    printf("  Total E0 = U_i =%.6e J\n", U_total);
   
    
    return U_total;
}

int main() {
    // Initialize 4 particles at square corners
    Particle particles[4];
    
    // Particle 1: Positron at (d/2, d/2)
    particles[0].x = d/2.0;
    particles[0].y = d/2.0;
    particles[0].vx = 0.0;
    particles[0].vy = 0.0;
    particles[0].mass = m_e;
    particles[0].charge = e;
    
    // Particle 2: Proton at (-d/2, d/2)
    particles[1].x = -d/2.0;
    particles[1].y = d/2.0;
    particles[1].vx = 0.0;
    particles[1].vy = 0.0;
    particles[1].mass = m_p;
    particles[1].charge = e;
    
    // Particle 3: Proton at (d/2, -d/2)
    particles[2].x = d/2.0;
    particles[2].y = -d/2.0;
    particles[2].vx = 0.0;
    particles[2].vy = 0.0;
    particles[2].mass = m_p;
    particles[2].charge = e;
    
    // Particle 4: Positron at (-d/2, -d/2)
    particles[3].x = -d/2.0;
    particles[3].y = -d/2.0;
    particles[3].vx = 0.0;
    particles[3].vy = 0.0;
    particles[3].mass = m_e;
    particles[3].charge = e;
    
    // Calculate initial energy
    double E0_calculated = calc_initial_potential_energy();
    double E0 = calculate_potential_energy(particles, 4);  
    printf("Verification: E0 from actual positions = %.6e J\n\n", E0);
    
    // Time parameters
    double dt = 1e-21;        
    double t = 0.0;
    double t_max = 1e-14;   
    int step = 0;
    int print_interval = 1000000; 
    
    printf("%-15s %-15s %-15s %-15s %-15s\n", "Time(s)", "KE(J)", "PE(J)", "Total_E(J)", "dE/E0");
    printf("---------------------------------------------------------------------------------\n");
    
    
    double KE_initial = calculate_kinetic_energy(particles, 4);
    double PE_initial = calculate_potential_energy(particles, 4);
    double E_initial = KE_initial + PE_initial;
    printf("%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n", 
           0.0, KE_initial, PE_initial, E_initial, 0.0);
    
    
    while (t < t_max) {
        
        rk2_step(particles, 4, dt);
        
        
        t += dt;
        step++;
        
        // Calculate energies
        double KE = calculate_kinetic_energy(particles, 4);
        double PE = calculate_potential_energy(particles, 4);
        double E_total = KE + PE;
        double relative_error = (E_total - E0) / E0;
        
        
        if (step % print_interval == 0) {
            printf("%-15.6e %-15.6e %-15.6e %-15.6e %-15.6e\n", 
                   t, KE, PE, E_total, relative_error);
        }
    }
    
  
    printf("\n=== Final State ===\n");
    printf("Time: %.6e s\n", t);
    
 
    double v_positron_total = 0.0, v_proton_total = 0.0;
    double KE_positron = 0.0, KE_proton = 0.0;
    int n_positron = 0, n_proton = 0;
    
    for (int i = 0; i < 4; i++) {
        double v = sqrt(particles[i].vx * particles[i].vx + particles[i].vy * particles[i].vy);
        double r = sqrt(particles[i].x * particles[i].x + particles[i].y * particles[i].y);
        double KE_i = 0.5 * particles[i].mass * v * v;
        
        printf("Particle %d: r=%.6e m, v=%.6e m/s, KE=%.6e J\n", i, r, v, KE_i);
        
        if (particles[i].mass == m_e) { 
            v_positron_total += v;
            KE_positron += KE_i;
            n_positron++;
        } else {  
            v_proton_total += v;
            KE_proton += KE_i;
            n_proton++;
        }
    }
    
    double v_positron_avg = v_positron_total / n_positron;
    double v_proton_avg = v_proton_total / n_proton;
    
    printf("\nSpeed ratio (positron/proton): %.6f\n", v_positron_avg / v_proton_avg);
    printf("KE ratio (positron/proton): %.6f\n", KE_positron / KE_proton);
    
    // Temperature calculation (kT = <KE> per particle)
    double T_positron = KE_positron / n_positron / (1.380649e-23);  // k_B
    double T_proton = KE_proton / n_proton / (1.380649e-23);
    printf("\n'Temperature' of positrons: %.6e K\n", T_positron);
    printf("'Temperature' of protons: %.6e K\n", T_proton);
    printf("Temperature ratio (positron/proton): %.6f\n", T_positron / T_proton);
    
    return 0;
}