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


// Function to estimate the fourth derivative using a higher-order finite difference formula
double fourth_derivative(double x[], double y[], int n, double* error) {
    if (n < 7) {
        printf("Not enough points to estimate the fourth derivative with higher accuracy.\n");
        *error = 0.0;
        return 0.0;
    }

    double max_d4 = 0.0;
    double max_error = 0.0;
    double h = x[1] - x[0]; // Assuming uniform spacing

    for (int i = 3; i < n - 3; i++) {
        // Higher-order finite difference formula for the fourth derivative
        double d4 = (-y[i-3] + 12 * y[i-2] - 39 * y[i-1] + 56 * y[i] - 39 * y[i+1] + 12 * y[i+2] - y[i+3]) / pow(h, 4);

        // Error estimation: proportional to the sixth derivative (unknown, so we approximate)
        double d4_error = fabs(d4) * pow(h, 2); // Error is O(h^2)

        // Track the maximum derivative and its associated error
        if (fabs(d4) > max_d4) {
            max_d4 = fabs(d4);
        }
        if (d4_error > max_error) {
            max_error = d4_error;
        }
    }

    *error = max_error;
    return max_d4;
}

/// Function to estimate the error of Simpson's rule
double simpson_error(double x[], double y[], int n) {
    if (n % 2 == 0) {
        printf("Number of intervals must be odd.\n");
        return 0.0;
    }

    double h = (x[n-1] - x[0]) / (n - 1); // Step size
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

    printf("Truncation Error: %e, Propagated Error: %e, Total Error: %e\n", truncation_error, propagated_error, total_error);

    return total_error;
}