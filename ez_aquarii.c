#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Constants
#define G 6.67430e-11  // Gravitational constant in m^3 kg^-1 s^-2
#define AU 1.496e11    // Astronomical Unit in meters
#define SOLAR_MASS 1.989e30  // Solar mass in kg
#define H 1000 // Time step size in seconds 
#define AVG_R_B 1.22 // Average separation of B from AC COM
#define E_B 0.439 // Eccentricity of B star orbit
#define NSTEPS 1000000 
#define OUTPUT_FREQ 1000 

// Global variables
int N_p;          // Number of particles
double *masses;   // Masses of particles

// Function prototypes
void deriv(double t, double y[], double dydx[]);
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
         void (*derivs)(double, double[], double[]));

int main(int argc, char *argv[]) {

    // Parse command-line arguments
    char *method = argv[1];

    // Number of particles in EZ Aquarii system
    N_p = 3;

    // Allocate memory for masses and state vectors
    int N = N_p * 6; // Total number of state variables
    masses = (double *)malloc(N_p * sizeof(double));
    double *y = (double *)malloc(N * sizeof(double));
    double *yout = (double *)malloc(N * sizeof(double));
    double *dydx = (double *)malloc(N * sizeof(double));

    if (masses == NULL || y == NULL || yout == NULL || dydx == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        exit(1);
    }

    // Define masses (in kg)
    masses[0] = 0.1187 * SOLAR_MASS; // EZ Aquarii A
    masses[1] = 0.1145 * SOLAR_MASS; // EZ Aquarii B
    masses[2] = 0.0930 * SOLAR_MASS; // EZ Aquarii C

    // Orbital periods (in seconds)
    double T_AC = 3.786516 * 24 * 3600;  // Inner binary period
    double T_B = 822.6 * 24 * 3600;      // B orbiting AC

    // Compute semi-major axes using Kepler's Third Law
    double M_AC = masses[0] + masses[2];  // Total mass of inner binary

    double a_AC = cbrt( G * M_AC * T_AC * T_AC / (4 * M_PI * M_PI) );  // Inner binary separation
    double a_B = AVG_R_B * AU * (1+ E_B);  // B's initial separation from AC COM @ apoapsis 

    // Initial positions and velocities

    // Star A (inner binary component)
    y[0] = -masses[2] / M_AC * a_AC;  // x position
    y[1] = 0.0;                       // y position
    y[2] = 0.0;                       // z position

    // Star B (outer star)
    y[6] = a_B;  // x position
    y[7] = 0.0;  // y position
    y[8] = 0.0;  // z position

    // Star C (inner binary component)
    y[12] = masses[0] / M_AC * a_AC;  // x position
    y[13] = 0.0;                      // y position
    y[14] = 0.0;                      // z position

    // Compute initial velocities 
    double v_A = masses[2] / M_AC * (2 * M_PI * a_AC / T_AC);
    double v_B = sqrt((G * (M_AC + masses[1])) * ((2/a_B)-(1/(1.22*AU))));
    double v_C = masses[0] / M_AC * (2 * M_PI * a_AC / T_AC);

    y[3] = 0.0;         // vx
    y[4] = v_A;         // vy
    y[5] = 0.0;         // vz

    y[15] = 0.0;        // vx
    y[16] = -v_C;       // vy
    y[17] = 0.0;        // vz

    y[9] = 0.0;         // vx
    y[10] = v_B;        // vy
    y[11] = 0.0;        // vz

    // Open output file
    FILE *output_file = fopen("ez_aquarii_output.txt", "w");
    if (!output_file) {
        fprintf(stderr, "Error: Cannot open output file.\n");
        exit(1);
    }

    // Write header to output file
    fprintf(output_file, "# time");
    for (int i = 0; i < N_p; i++) {
        fprintf(output_file, " x%d y%d z%d vx%d vy%d vz%d", i + 1, i + 1, i + 1, i + 1, i + 1, i + 1);
    }
    fprintf(output_file, " placehold\n");

    // Main integration loop
    double t = 0.0;
    for (int step = 1; step <= NSTEPS; step++) {
        // Compute derivatives
        deriv(t, y, dydx);

        // Advance the system 
        if (strcmp(method, "rk4") == 0) {
            rk4(y, dydx, N, t, H, yout, deriv);
        } else {
            fprintf(stderr, "Error: Unknown method '%s'\n", method);
            exit(1);
        }

        // Update time and state
        t += H;
        for (int i = 0; i < N; i++) {
            y[i] = yout[i];
        }

        // Output results 
        if (step % OUTPUT_FREQ == 0) {

            // Output time, positions, velocities
            fprintf(output_file, "%e", t);
            for (int i = 0; i < N_p; i++) {
                fprintf(output_file, " %e %e %e %e %e %e",
                        y[i * 6 + 0], y[i * 6 + 1], y[i * 6 + 2],
                        y[i * 6 + 3], y[i * 6 + 4], y[i * 6 + 5]);
            }
            fprintf(output_file, " %e\n", 0.0000);
        }
    }

    // Clean up
    fclose(output_file);
    free(y);
    free(yout);
    free(dydx);
    free(masses);

    return 0;
}

// Computes derivatives of positions and velocities
void deriv(double t, double y[], double dydx[]) {

    // Initialize derivatives
    for (int i = 0; i < N_p; i++) {
        // Positions derivatives
        dydx[i * 6 + 0] = y[i * 6 + 3]; // dx/dt = vx
        dydx[i * 6 + 1] = y[i * 6 + 4]; // dy/dt = vy
        dydx[i * 6 + 2] = y[i * 6 + 5]; // dz/dt = vz

        // Initialize accelerations to zero
        double ax = 0.0, ay = 0.0, az = 0.0;

        // Compute gravitational acceleration
        for (int j = 0; j < N_p; j++) {
            if (j != i) {
                double dx = y[j * 6 + 0] - y[i * 6 + 0];
                double dy = y[j * 6 + 1] - y[i * 6 + 1];
                double dz = y[j * 6 + 2] - y[i * 6 + 2];
                double dist2 = dx * dx + dy * dy + dz * dz;
                double invDist3 = pow(dist2, -1.5);
                double force = G * masses[j] * invDist3;

                ax += force * dx;
                ay += force * dy;
                az += force * dz;
            }
        }

        // Velocities derivatives 
        dydx[i * 6 + 3] = ax;
        dydx[i * 6 + 4] = ay;
        dydx[i * 6 + 5] = az;
    }
}

// Runge-Kutta 4th order integration method using double 
void rk4(double y[], double dydx[], int n, double x, double h, double yout[],
         void (*derivs)(double, double[], double[])) {
    double *dym = (double *)malloc(n * sizeof(double));
    double *dyt = (double *)malloc(n * sizeof(double));
    double *yt = (double *)malloc(n * sizeof(double));

    double hh = h * 0.5;
    double h6 = h / 6.0;
    double xh = x + hh;

    // First step
    for (int i = 0; i < n; i++) {
        yt[i] = y[i] + hh * dydx[i];
    }

    // Second step
    (*derivs)(xh, yt, dyt);
    for (int i = 0; i < n; i++) {
        yt[i] = y[i] + hh * dyt[i];
    }

    // Third step
    (*derivs)(xh, yt, dym);
    for (int i = 0; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
        dym[i] += dyt[i];
    }

    // Fourth step
    (*derivs)(x + h, yt, dyt);
    for (int i = 0; i < n; i++) {
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    }

    free(dym);
    free(dyt);
    free(yt);
}

