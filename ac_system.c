// A-C system motion

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// Constants
#define SOLAR_MASS 1.989e30        // Solar mass in kg
#define AU 1.496e11                // Astronomical Unit in meters
#define G 6.67430e-11              // Gravitational constant in m^3 kg^-1 s^-2

// Star Structure
typedef struct {
    double m;    // Mass in kg
    double x;    // Position X in meters
    double y;    // Position Y in meters
    double z;    // Position Z in meters
    double vx;   // Velocity X in m/s
    double vy;   // Velocity Y in m/s
    double vz;   // Velocity Z in m/s
} Star;

// Function Prototypes
void compute_accelerations(Star *A, Star *B, Star *C,
                           double *ax_A, double *ay_A, double *az_A,
                           double *ax_B, double *ay_B, double *az_B,
                           double *ax_C, double *ay_C, double *az_C);
void runge_kutta4(Star *A, Star *B, Star *C, double dt);
void leapfrog_step(Star *A, Star *B, Star *C, double dt);

// Gravitational acceleration functions
void compute_accelerations(Star *A, Star *B, Star *C,
                           double *ax_A, double *ay_A, double *az_A,
                           double *ax_B, double *ay_B, double *az_B,
                           double *ax_C, double *ay_C, double *az_C) {
    // Reset accelerations
    *ax_A = *ay_A = *az_A = 0.0;
    *ax_C = *ay_C = *az_C = 0.0;

    // Acceleration on A due to C
    double dx_AC = C->x - A->x;
    double dy_AC = C->y - A->y;
    double dz_AC = C->z - A->z;
    double r_AC = sqrt(dx_AC*dx_AC + dy_AC*dy_AC + dz_AC*dz_AC);
    if (r_AC != 0) {
        double a_AC = G * C->m / (r_AC * r_AC * r_AC);
        *ax_A += a_AC * dx_AC;
        *ay_A += a_AC * dy_AC;
        *az_A += a_AC * dz_AC;
    }

    // Acceleration on C due to A
    double dx_CA = A->x - C->x;
    double dy_CA = A->y - C->y;
    double dz_CA = A->z - C->z;
    double r_CA = sqrt(dx_CA*dx_CA + dy_CA*dy_CA + dz_CA*dz_CA);
    if (r_CA != 0) {
        double a_CA = G * A->m / (r_CA * r_CA * r_CA);
        *ax_C += a_CA * dx_CA;
        *ay_C += a_CA * dy_CA;
        *az_C += a_CA * dz_CA;
    }
}


// Runge-Kutta 4 (RK4) Integration Method 
void runge_kutta4(Star *A, Star *B, Star *C, double dt) {
    // Store original states
    Star A_orig = *A;
    Star B_orig = *B;
    Star C_orig = *C;

    // Compute k1 accelerations
    double ax1_A, ay1_A, az1_A;
    double ax1_B, ay1_B, az1_B;
    double ax1_C, ay1_C, az1_C;
    compute_accelerations(A, B, C,
                          &ax1_A, &ay1_A, &az1_A,
                          &ax1_B, &ay1_B, &az1_B,
                          &ax1_C, &ay1_C, &az1_C);

    // Compute k1
    Star k1_A = {
        A->m,
        A->x,
        A->y,
        A->z,
        A->vx,
        A->vy,
        A->vz
    };
    Star k1_B = {
        B->m,
        B->x,
        B->y,
        B->z,
        B->vx,
        B->vy,
        B->vz
    };
    Star k1_C = {
        C->m,
        C->x,
        C->y,
        C->z,
        C->vx,
        C->vy,
        C->vz
    };

    // Estimate mid-step stars for k2
    Star A_mid = {
        A->m,
        A->x + 0.5 * dt * A->vx,
        A->y + 0.5 * dt * A->vy,
        A->z + 0.5 * dt * A->vz,
        A->vx + 0.5 * dt * ax1_A,
        A->vy + 0.5 * dt * ay1_A,
        A->vz + 0.5 * dt * az1_A
    };
    Star B_mid = {
        B->m,
        B->x + 0.5 * dt * B->vx,
        B->y + 0.5 * dt * B->vy,
        B->z + 0.5 * dt * B->vz,
        B->vx + 0.5 * dt * ax1_B,
        B->vy + 0.5 * dt * ay1_B,
        B->vz + 0.5 * dt * az1_B
    };
    Star C_mid = {
        C->m,
        C->x + 0.5 * dt * C->vx,
        C->y + 0.5 * dt * C->vy,
        C->z + 0.5 * dt * C->vz,
        C->vx + 0.5 * dt * ax1_C,
        C->vy + 0.5 * dt * ay1_C,
        C->vz + 0.5 * dt * az1_C
    };

    // Compute k2 accelerations
    double ax2_A, ay2_A, az2_A;
    double ax2_B, ay2_B, az2_B;
    double ax2_C, ay2_C, az2_C;
    compute_accelerations(&A_mid, &B_mid, &C_mid,
                          &ax2_A, &ay2_A, &az2_A,
                          &ax2_B, &ay2_B, &az2_B,
                          &ax2_C, &ay2_C, &az2_C);

    // Compute k2
    Star k2_A = {
        A_mid.m,
        A_mid.x,
        A_mid.y,
        A_mid.z,
        A_mid.vx,
        A_mid.vy,
        A_mid.vz
    };
    Star k2_B = {
        B_mid.m,
        B_mid.x,
        B_mid.y,
        B_mid.z,
        B_mid.vx,
        B_mid.vy,
        B_mid.vz
    };
    Star k2_C = {
        C_mid.m,
        C_mid.x,
        C_mid.y,
        C_mid.z,
        C_mid.vx,
        C_mid.vy,
        C_mid.vz
    };

    // Estimate mid-step stars for k3
    Star A_mid2 = {
        A->m,
        A->x + 0.5 * dt * A_mid.vx,
        A->y + 0.5 * dt * A_mid.vy,
        A->z + 0.5 * dt * A_mid.vz,
        A->vx + 0.5 * dt * ax2_A,
        A->vy + 0.5 * dt * ay2_A,
        A->vz + 0.5 * dt * az2_A
    };
    Star B_mid2 = {
        B->m,
        B->x + 0.5 * dt * B_mid.vx,
        B->y + 0.5 * dt * B_mid.vy,
        B->z + 0.5 * dt * B_mid.vz,
        B->vx + 0.5 * dt * ax2_B,
        B->vy + 0.5 * dt * ay2_B,
        B->vz + 0.5 * dt * az2_B
    };
    Star C_mid2 = {
        C->m,
        C->x + 0.5 * dt * C_mid.vx,
        C->y + 0.5 * dt * C_mid.vy,
        C->z + 0.5 * dt * C_mid.vz,
        C->vx + 0.5 * dt * ax2_C,
        C->vy + 0.5 * dt * ay2_C,
        C->vz + 0.5 * dt * az2_C
    };

    // Compute k3 accelerations
    double ax3_A, ay3_A, az3_A;
    double ax3_B, ay3_B, az3_B;
    double ax3_C, ay3_C, az3_C;
    compute_accelerations(&A_mid2, &B_mid2, &C_mid2,
                          &ax3_A, &ay3_A, &az3_A,
                          &ax3_B, &ay3_B, &az3_B,
                          &ax3_C, &ay3_C, &az3_C);

    // Compute k3
    Star k3_A = {
        A_mid2.m,
        A_mid2.x,
        A_mid2.y,
        A_mid2.z,
        A_mid2.vx,
        A_mid2.vy,
        A_mid2.vz
    };
    Star k3_B = {
        B_mid2.m,
        B_mid2.x,
        B_mid2.y,
        B_mid2.z,
        B_mid2.vx,
        B_mid2.vy,
        B_mid2.vz
    };
    Star k3_C = {
        C_mid2.m,
        C_mid2.x,
        C_mid2.y,
        C_mid2.z,
        C_mid2.vx,
        C_mid2.vy,
        C_mid2.vz
    };

    // Estimate end-step stars for k4
    Star A_end = {
        A->m,
        A->x + dt * A_mid.vx,
        A->y + dt * A_mid.vy,
        A->z + dt * A_mid.vz,
        A->vx + dt * ax3_A,
        A->vy + dt * ay3_A,
        A->vz + dt * az3_A
    };
    Star B_end = {
        B->m,
        B->x + dt * B_mid.vx,
        B->y + dt * B_mid.vy,
        B->z + dt * B_mid.vz,
        B->vx + dt * ax3_B,
        B->vy + dt * ay3_B,
        B->vz + dt * az3_B
    };
    Star C_end = {
        C->m,
        C->x + dt * C_mid.vx,
        C->y + dt * C_mid.vy,
        C->z + dt * C_mid.vz,
        C->vx + dt * ax3_C,
        C->vy + dt * ay3_C,
        C->vz + dt * az3_C
    };

    // Compute k4 accelerations
    double ax4_A, ay4_A, az4_A;
    double ax4_B, ay4_B, az4_B;
    double ax4_C, ay4_C, az4_C;
    compute_accelerations(&A_end, &B_end, &C_end,
                          &ax4_A, &ay4_A, &az4_A,
                          &ax4_B, &ay4_B, &az4_B,
                          &ax4_C, &ay4_C, &az4_C);

    // Compute k4
    Star k4_A = {
        A_end.m,
        A_end.x,
        A_end.y,
        A_end.z,
        A_end.vx,
        A_end.vy,
        A_end.vz
    };
    Star k4_B = {
        B_end.m,
        B_end.x,
        B_end.y,
        B_end.z,
        B_end.vx,
        B_end.vy,
        B_end.vz
    };
    Star k4_C = {
        C_end.m,
        C_end.x,
        C_end.y,
        C_end.z,
        C_end.vx,
        C_end.vy,
        C_end.vz
    };

    // Update positions and velocities using RK4 formula
    A->x = A_orig.x + (dt / 6.0) * (A_mid.vx + 2.0 * A_mid2.vx + A_end.vx);
    A->y = A_orig.y + (dt / 6.0) * (A_mid.vy + 2.0 * A_mid2.vy + A_end.vy);
    A->z = A_orig.z + (dt / 6.0) * (A_mid.vz + 2.0 * A_mid2.vz + A_end.vz);
    A->vx = A_orig.vx + (dt / 6.0) * (ax1_A + 2.0 * ax2_A + 2.0 * ax3_A + ax4_A);
    A->vy = A_orig.vy + (dt / 6.0) * (ay1_A + 2.0 * ay2_A + 2.0 * ay3_A + ay4_A);
    A->vz = A_orig.vz + (dt / 6.0) * (az1_A + 2.0 * az2_A + 2.0 * az3_A + az4_A);

    B->x = B_orig.x + (dt / 6.0) * (B_mid.vx + 2.0 * B_mid2.vx + B_end.vx);
    B->y = B_orig.y + (dt / 6.0) * (B_mid.vy + 2.0 * B_mid2.vy + B_end.vy);
    B->z = B_orig.z + (dt / 6.0) * (B_mid.vz + 2.0 * B_mid2.vz + B_end.vz);
    B->vx = B_orig.vx + (dt / 6.0) * (ax1_B + 2.0 * ax2_B + 2.0 * ax3_B + ax4_B);
    B->vy = B_orig.vy + (dt / 6.0) * (ay1_B + 2.0 * ay2_B + 2.0 * ay3_B + ay4_B);
    B->vz = B_orig.vz + (dt / 6.0) * (az1_B + 2.0 * az2_B + 2.0 * az3_B + az4_B);

    C->x = C_orig.x + (dt / 6.0) * (C_mid.vx + 2.0 * C_mid2.vx + C_end.vx);
    C->y = C_orig.y + (dt / 6.0) * (C_mid.vy + 2.0 * C_mid2.vy + C_end.vy);
    C->z = C_orig.z + (dt / 6.0) * (C_mid.vz + 2.0 * C_mid2.vz + C_end.vz);
    C->vx = C_orig.vx + (dt / 6.0) * (ax1_C + 2.0 * ax2_C + 2.0 * ax3_C + ax4_C);
    C->vy = C_orig.vy + (dt / 6.0) * (ay1_C + 2.0 * ay2_C + 2.0 * ay3_C + ay4_C);
    C->vz = C_orig.vz + (dt / 6.0) * (az1_C + 2.0 * az2_C + 2.0 * az3_C + az4_C);
}

// Leapfrog Integration Method
void leapfrog_step(Star *A, Star *B, Star *C, double dt) {
    // Compute initial accelerations
    double ax_A, ay_A, az_A;
    double ax_B, ay_B, az_B;
    double ax_C, ay_C, az_C;
    compute_accelerations(A, B, C,
                          &ax_A, &ay_A, &az_A,
                          &ax_B, &ay_B, &az_B,
                          &ax_C, &ay_C, &az_C);

    // Half-step velocity update
    A->vx += 0.5 * dt * ax_A;
    A->vy += 0.5 * dt * ay_A;
    A->vz += 0.5 * dt * az_A;

    B->vx += 0.5 * dt * ax_B;
    B->vy += 0.5 * dt * ay_B;
    B->vz += 0.5 * dt * az_B;

    C->vx += 0.5 * dt * ax_C;
    C->vy += 0.5 * dt * ay_C;
    C->vz += 0.5 * dt * az_C;

    // Full-step position update
    A->x += dt * A->vx;
    A->y += dt * A->vy;
    A->z += dt * A->vz;

    B->x += dt * B->vx;
    B->y += dt * B->vy;
    B->z += dt * B->vz;

    C->x += dt * C->vx;
    C->y += dt * C->vy;
    C->z += dt * C->vz;

    // Compute new accelerations after position update
    compute_accelerations(A, B, C,
                          &ax_A, &ay_A, &az_A,
                          &ax_B, &ay_B, &az_B,
                          &ax_C, &ay_C, &az_C);

    // Half-step velocity update with new accelerations
    A->vx += 0.5 * dt * ax_A;
    A->vy += 0.5 * dt * ay_A;
    A->vz += 0.5 * dt * az_A;

    B->vx += 0.5 * dt * ax_B;
    B->vy += 0.5 * dt * ay_B;
    B->vz += 0.5 * dt * az_B;

    C->vx += 0.5 * dt * ax_C;
    C->vy += 0.5 * dt * ay_C;
    C->vz += 0.5 * dt * az_C;
}

int main(){
    // Convert solar masses to kg
    double mass_A = 0.1187 * SOLAR_MASS; // EZ Aquarii A
    double mass_B = 0.1145 * SOLAR_MASS; // EZ Aquarii B
    double mass_C = 0.0930 * SOLAR_MASS; // EZ Aquarii C

    // Orbital parameters
    double separation_AC = 0.03 * AU; // Total separation of A and C in meters
    double separation_AC_A = (mass_C / (mass_A + mass_C)) * separation_AC; // Separation of A from COM
    double separation_AC_C = (mass_A / (mass_A + mass_C)) * separation_AC; // Separation of C from COM
    double separation_B_AC = 1.5 * AU; // Approximate initial separation of B from AC pair

    // Initial positions 
    Star A = {mass_A, -separation_AC_A, 0.0, 0.0, 0.0, 0.0, 0.0};
    Star C = {mass_C, separation_AC_C, 0.0, 0.0, 0.0, 0.0, 0.0};
    Star B = {mass_B, 0.0, separation_B_AC, 0.0, 0.0, 0.0, 0.0};

    // Calculate initial velocities based on orbital periods
    // Convert periods from days to seconds
    double T_AC = 3.786516 * 24 * 3600; // 3.786516 days in seconds
    double T_B_AC = 823.6 * 24 * 3600;  // 823.6 days in seconds

    // Velocities for A and C to orbit each other (circular orbit)
    double v_A_initial = 2.0 * M_PI * separation_AC_A / T_AC;
    double v_C_initial = 2.0 * M_PI * separation_AC_C / T_AC;

    // Assign velocities 
    A.vy = v_A_initial;    // Star A velocity (positive y-direction)
    C.vy = -v_C_initial;   // Star C velocity (negative y-direction)

    // Velocity for B to orbit the AC pair
    double v_B_initial = 2.0 * M_PI * separation_B_AC / T_B_AC;
    B.vx = v_B_initial;    // Star B velocity (positive x-direction)

    // Open output files for RK4 method
    FILE *rk4_file = fopen("ac_output.txt", "w");

    if (rk4_file == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    // Write headers 
    fprintf(rk4_file, "Timestep\tA_X\tA_Y\tA_Z\tA_VX\tA_VY\tA_VZ\tB_X\tB_Y\tB_Z\tB_VX\tB_VY\tB_VZ\tC_X\tC_Y\tC_Z\tC_VX\tC_VY\tC_VZ\n");

    // Time step and number of simulation steps
    double dt = 60; // Time step in seconds (1 minute)
    int steps = 10000; // Total simulation time

    // Simulation loop
    for (int i = 0; i < steps; i++) { 
        double current_time = dt * (i + 1); // Current simulation time in seconds

        leapfrog_step(&A, &B, &C, dt);

        // Runge-Kutta 4 Integration
        runge_kutta4(&A, &B, &C, dt);

        // Write RK4 results to file
        fprintf(rk4_file, "%.2f\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\n", 
                current_time, 
                A.x, A.y, A.z, 
                A.vx, A.vy, A.vz, 
                B.x, B.y, B.z, 
                B.vx, B.vy, B.vz,
                C.x, C.y, C.z,
                C.vx, C.vy, C.vz);
    }

    // Close output files
    fclose(rk4_file);

    return 0;
}