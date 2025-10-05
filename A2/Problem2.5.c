#include <stdio.h>
#include <math.h>

// Function to estimate distance using Titius-Blagg law: d = x^n
double estimate_distance(int n, double x){
    return pow(x, n);
}

// Function to calculate the quadratic sum of relative errors: E = sum((d_i - a_i)/a_i)^2
double calculate_error(double x, double a[], int n_values[], int size){
    double total_error = 0.0;
    
    for(int i = 0; i < size; i++){
        double predicted = estimate_distance(n_values[i], x);
        double relative_error = pow((predicted - a[i]) / a[i],2);
        total_error += relative_error;
    }
    
    return total_error;
}

// Derivative of E with respect to x (for bisection method)
double error_derivative(double x, double a[], int n[], int size){
    /* we first find the derivate on paper and then use the function we get to estimate the derivative of the relative error*/
    double derivative = 0.0;
    
    for(int i = 0; i < size; i++){
        double dE_dx= 2*((estimate_distance(n[i], x) - a[i]) / a[i]) * n[i] * pow(x, n[i] - 1);
        derivative += dE_dx;
    }
    
    return derivative;
}

// Bisection method to find optimal x
double find_optimal_x(double a[], int n_values[], int size, double x_min, double x_max, double tolerance){
    double x_left = x_min;
    double x_right = x_max;
    
    printf("Finding optimal x using bisection method:\n");
    printf("Initial range: [%.3f, %.3f]\n\n", x_left, x_right);
    // we assume that the bisection method will converge in 50 iterations
    for(int iteration = 0; iteration < 50; iteration++){
        double x_mid = (x_left + x_right) / 2.0;
        double derivative_mid = error_derivative(x_mid, a, n_values, size);
        double derivative_left = error_derivative(x_left, a, n_values, size);
        double derivative_right = error_derivative(x_right, a, n_values, size);
        
        printf("Iteration %2d: x_mid = %.6f, derivative = %.6f\n", 
               iteration + 1, x_mid, derivative_mid);
        
        if(fabs(derivative_mid) < tolerance){
            printf("\nConverged! Optimal x = %.6f\n", x_mid);
            return x_mid;
        }
        
        if(derivative_mid * derivative_left < 0){
            x_right = x_mid;
        } else {
            x_left = x_mid;
        }
    }
    
    return (x_left + x_right) / 2.0;
}

int main(){
    // Semi-major axis distances in AU (observational data)-> true value of the distances of planets fro sun 
    double a[9] = {
        0.3871,      // Mercury (n = -2)
        0.7233,      // Venus   (n = -1)  
        1.0000,      // Earth   (n = 0)
        1.5236,      // Mars    (n = 1)
        2.7692,      // Ceres   (n = 2)
        5.2034,      // Jupiter (n = 3)
        9.5371,      // Saturn  (n = 4)
        19.165,      // Uranus  (n = 5)
        30.178       // Neptune (n = 6)
    };
    
    // The n values from -2 to 6 for the 9 celestial bodies
    int n[9] = {-2, -1, 0, 1, 2, 3, 4, 5, 6};
    
    // Names for output
    char* names[9] = {
        "Mercury", "Venus", "Earth", "Mars", "Ceres", 
        "Jupiter", "Saturn", "Uranus", "Neptune"
    };
 
    
    printf("Observational data (semi-major axis in AU):\n");
    for(int i = 0; i < 9; i++){
        printf("%-8s (n=%2d): %7.4f AU\n", names[i], n[i], a[i]);
    }
    printf("\n");
    
    // Find optimal x using bisection method
    double optimal_x = find_optimal_x(a, n, 9, 1.5, 2.5, 1e-6);
    
    printf("\n" "Results with optimal x = %.6f:\n", optimal_x);
   
    // Use our existing function to calculate total error
    double total_error = calculate_error(optimal_x, a, n, 9);
    
    for(int i = 0; i < 9; i++){
        double predicted = estimate_distance(n[i], optimal_x);
        double relative_error = (predicted - a[i]) / a[i];
        double absolute_error = predicted - a[i];
        
        printf(" The relative_error, absolute error and the predicted distance for the optimal x value for all planets is: %-8s %2d   %7.4f    %7.4f    %8.3f%%  %8.4f\n", 
               names[i], n[i], a[i], predicted, 
               relative_error * 100, absolute_error);
    }
    
    printf("\nTotal quadratic error E = %.6f\n", total_error);
    printf("\n The optimal x value for the Titius-Blagg law: distance = %.6f^n AU\n", optimal_x);
    
    return 0;
}
    

   
