#include <stdio.h>
#include <math.h>

// Function to calculate the integral using Simpson's rule
double simpson(double x[], double y[], int n) {
    if (n % 2 == 0) {
        printf("Number of intervals must be odd.\n");
        return 0.0;
    }

    double h = (x[n-1] - x[0]) / (n - 1);
    double sum = y[0] + y[n-1];

    for (int i = 1; i < n-1; i += 2) {
        sum += 4 * y[i];
    }

    for (int i = 2; i < n-2; i += 2) {
        sum += 2 * y[i];
    }

    return (h / 3) * sum;
}


// Function to estimate the fourth derivative using finite differences
double fourth_derivative(double x[], double y[], int n) {
    if (n < 5) {
        printf("Not enough points to estimate the fourth derivative.\n");
        return 0.0;
    }

    double max_d4 = 0.0;
    for (int i = 2; i < n - 2; i++) {
        double d4 = (y[i-2] - 4 * y[i-1] + 6 * y[i] - 4 * y[i+1] + y[i+2]) / pow(x[1] - x[0], 4);
        if (fabs(d4) > max_d4) {
            max_d4 = fabs(d4);
        }
    }
    return max_d4;
}

// Function to estimate the error of Simpson's rule
double simpson_error(double x[], double y[], int n) {
    if (n % 2 == 0) {
        printf("Number of intervals must be odd.\n");
        return 0.0;
    }

    double h = (x[n-1] - x[0]) / (n - 1);
    double max_d4 = fourth_derivative(x, y, n);

    // Error formula: E = -(b-a)^5 / (180 * n^4) * max|f^(4)(xi)|
    double error = (pow(x[n-1] - x[0], 5) / (180 * pow(n - 1, 4))) * max_d4;
    return error;
}