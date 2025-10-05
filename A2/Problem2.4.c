#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <time.h>

double max(double t1, double t2){
    if(t1 > t2){
        return t1;
    }
    else{
        return t2;
    }
}

double min(double t1, double t2){
    if(t1 < t2){
        return t1;
    }
    else{
        return t2;
    }
}

double uniform(){
    return (double)rand() / (double)RAND_MAX;
}

double absolute_value(double t){
    if(t < 0){
        return -t;
    }
    else{
        return t;
    }
}

int main(){
    srand(time(NULL));
    int N[]= {100, 1000, 10000, 100000, 1000000,10000000,100000000};
    int length_N= sizeof(N) / sizeof(N[0]);
    
    // True values for the probability, minimum, absolute value and maximum to calculate the relative errors
    double TRUE_PROBABILITY = 4.0 / 9.0;
    double TRUE_MIN = 1.0 / 3.0;
    double TRUE_ABS = 1.0 / 3.0;
    double TRUE_MAX = 1.0 / 3.0;
    
    // Open file for convergence data
    FILE *fp = fopen("convergence.dat", "w");
    fprintf(fp, "# N    Rel_Error_Prob    Rel_Error_Min    Rel_Error_Abs    Rel_Error_Max\n");
    
    for (int i = 0; i < length_N; i++){
        //in the first part we had to calculate the probablity of P(|t1-t2|>1/3)
        // to calculate that we first need to sample from the uniform distribution N times and find the ratio
        // of the number of times that the absolute value of the difference between t1 and t2 is greater than 1/3 to N.
        // In the second part we need to calculate the following average time spans: (i) before the arrival of the first particle, (ii) between their arrivals, (iii) after the second arrival
        // These expected values are given by (i) min(t1,t2), (ii)E(|t1-t2|), (iii) E(1-max(t1,t2))

        int count= 0;
        double sum_min= 0;
        double sum_absolute_value= 0;
        double sum_max= 0;
        // We also have to calculate the standard deviation of the three expected values
        // r^2= <r^2> - <r>^2
        double min_value_squared=0;
        double absolute_value_value_squared=0;
        double max_value_squared=0;

        for (int j = 0; j < N[i]; j++){
            double t1= uniform();
            double t2= uniform();
            if (absolute_value(t1-t2)> 1.0/3.0){
                count++;
            }
            sum_min+= min(t1,t2);
            sum_absolute_value+= absolute_value(t1-t2);
            sum_max+= max(t1,t2);
            min_value_squared+= min(t1,t2)*min(t1,t2);
            absolute_value_value_squared+= absolute_value(t1-t2)*absolute_value(t1-t2);
            max_value_squared+= max(t1,t2)*max(t1,t2);
        }
        double probability= (double)count / N[i];
        double expected_value_min= sum_min / N[i];
        double expected_value_absolute_value= sum_absolute_value / N[i];
        double expected_value_max= sum_max / N[i];
        double expected_value_min_squared= min_value_squared / N[i];
        double expected_value_absolute_value_squared= absolute_value_value_squared / N[i];
        double expected_value_max_squared= max_value_squared / N[i];
        double standard_deviation_min= sqrt(expected_value_min_squared - expected_value_min*expected_value_min);
        double standard_deviation_absolute_value= sqrt(expected_value_absolute_value_squared - expected_value_absolute_value*expected_value_absolute_value);
        double standard_deviation_max= sqrt(expected_value_max_squared - expected_value_max*expected_value_max);
        
        // Calculate relative errors for convergence analysis
        double rel_error_prob = fabs(probability - TRUE_PROBABILITY) / TRUE_PROBABILITY;
        double rel_error_min = fabs(expected_value_min - TRUE_MIN) / TRUE_MIN;
        double rel_error_abs = fabs(expected_value_absolute_value - TRUE_ABS) / TRUE_ABS;
        double rel_error_max = fabs((1.0 - expected_value_max) - TRUE_MAX) / TRUE_MAX;
        
        // This line writes the relative errors to the convergence.dat file
        fprintf(fp, "%d    %.10e    %.10e    %.10e    %.10e\n", 
                N[i], rel_error_prob, rel_error_min, rel_error_abs, rel_error_max);
        
        printf("/****************************/\n");
        printf("The probability of P(|t1-t2|>1/3) is %f for N=%d\n", probability, N[i]);
        printf("The expected value of min(t1,t2) is %f for N=%d\n", expected_value_min, N[i]);
        printf("The expected value of E(|t1-t2|) is %f for N=%d\n", expected_value_absolute_value, N[i]);
        printf("The expected value of E(1-max(t1,t2)) is %f for N=%d\n", 1-expected_value_max, N[i]);
        printf("The standard deviation of min(t1,t2) is %f for N=%d\n", standard_deviation_min, N[i]);
        printf("The standard deviation of E(|t1-t2|) is %f for N=%d\n", standard_deviation_absolute_value, N[i]);
        printf("The standard deviation of E(1-max(t1,t2)) is %f for N=%d\n", standard_deviation_max, N[i]);
        printf("Relative error (probability): %.6e\n", rel_error_prob);
        printf("/****************************/\n");
        printf("\n");
    }
    
    fclose(fp);
   
    
    return 0;

}