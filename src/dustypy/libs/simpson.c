#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct {
    double value;
    double error;
} SimpsonResult;


// Function to estimate the fourth derivative using a higher-order finite difference formula
double fourth_derivative(double x[], double y[], int n, double* error_der) {
    if (n < 7) { *error_der = 0.0; return 0.0; }
    double max_d4 = 0.0;
    double h = x[1] - x[0];
    double h4 = pow(h, 4);

    for (int i = 3; i < n - 3; i++) {
        double d4 = (-y[i-3] + 12*y[i-2] - 39*y[i-1] + 56*y[i] - 39*y[i+1] + 12*y[i+2] - y[i+3]) / h4;
        if (fabs(d4) > max_d4) max_d4 = fabs(d4);
    }
    *error_der = max_d4 * pow(h, 2);
    return max_d4;
}

/// Function to estimate the error of Simpson's rule
double simpson_error(double x[], double y[], int n) {
    if (n % 2 == 0) {
        printf("Number of intervals must be odd.\n");
        return 0.0;
    }

    double a = x[0];                     // Start of the interval
    double b = x[n-1];                   // End of the interval
    double error_der;

    // Calculate the maximum fourth derivative and its error
    double max_d4 = fourth_derivative(x, y, n, &error_der);

    // Error formula: E = -(b-a)^5 / (180 * n^4) * max|f^(4)(xi)|
    double truncation_error = fabs((pow(b - a, 5) / (180 * pow(n - 1, 4))) * max_d4);

    // Propagate the error of the fourth derivative into the Simpson error
    double propagated_error = fabs((pow(b - a, 5) / (180 * pow(n - 1, 4))) * error_der);

    // Combine the truncation error and propagated error
    double total_error = sqrt(pow(truncation_error, 2) + pow(propagated_error, 2));
    
    return total_error;
}

// Function to calculate the integral using Simpson's rule
SimpsonResult simpson(double x[], double y[], int n) {
    SimpsonResult res = {0.0, 0.0};

    if (n % 2 == 0) {
        printf("Number of intervals must be odd.\n");
        return res;
    }

    double y_norm[n];

    double max_y = 0.0;
    for (int i = 0; i < n; i++) {
        if (fabs(y[i]) > max_y) max_y = fabs(y[i]);
    }

    for (int i = 0; i < n; i++) {
        y_norm[i] = y[i]/max_y;
    }

    double h = (x[n-1] - x[0]) / (n - 1);
    double sum = y_norm[0] + y_norm[n-1];

    #pragma omp parallel for reduction(+:sum)
    for (int i = 1; i < n-1; i += 2) {
        sum += 4 * y_norm[i];
    }

    #pragma omp parallel for reduction(+:sum)
    for (int i = 2; i < n-2; i += 2) {
        sum += 2 * y_norm[i];
    }

    res.value = (h / 3) * sum * max_y;
    res.error = simpson_error(x, y_norm, n) * max_y;

    return res;
}
