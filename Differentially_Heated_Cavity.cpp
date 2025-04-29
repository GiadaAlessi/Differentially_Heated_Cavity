#include "Differentially_Heated_Cavity.h"
#include<iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ranges>

using namespace std;
using namespace std::chrono;

// Function to allocate matrices
vector<vector<double>> CreateMatrix(int rows, int cols, double init_val = 0.0) {
    return vector(rows, vector(cols, init_val));
}

// Base Class to implement multiple methods
class ConvectiveScheme {
public:
    virtual void ComputeRt(double rho, int Ny, int Nx, double dx, double lambda, double cp, vector<vector<double> > &u,
                           vector<vector<double> > &v, vector<vector<double> > &Rt,
                           vector<vector<double> > &T) = 0; // Pure virtual function

    virtual void ComputeRu(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Ru,
                           double rho, double dx, double mu, int Nx, int Ny) = 0;

    virtual void ComputeRv(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Rv,
                           double rho, double dx, double mu, int Nx, int Ny) = 0;

    virtual void Nusselt(int Nx, int Ny, double Tcold, double Thot, double L, double alpha, double dx, double &Nux_avg,
                         vector<vector<double> > &T1star, vector<vector<double> > &T1, vector<vector<double> > &qx,
                         vector<vector<double> > &u1star, vector<vector<double> > &u1, vector<vector<double> > &dTdx,
                         vector<double> &Nux, double dxstar) = 0;

    virtual ~ConvectiveScheme() = default; // Virtual destructor
};


// Central Differencing Scheme (CDS) method
class CDSMethod : public ConvectiveScheme {

public:
    // Function to compute R(t) on each node of stagg-T mesh (same as stagg-P)
    void ComputeRt(double rho, int Ny, int Nx, double dx, double lambda, double cp, vector<vector<double> > &u,
                   vector<vector<double> > &v, vector<vector<double> > &Rt, vector<vector<double> > &T) override {
        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                Rt[i][j] = -rho * dx * cp * (u[i][j + 1] * 1.0 / 2 * (T[i][j] + T[i][j + 1]) - u[i][j] * 1.0 / 2 * (
                                            T[i][j] + T[i][j - 1]) + v[i][j] * 1.0 / 2 * (T[i][j] + T[i - 1][j]) - v[
                                            i + 1][j] * 1.0 / 2 * (T[i][j] + T[i + 1][j])) + lambda * (
                               T[i][j + 1] + T[i][j - 1] + T[i - 1][j] + T[i + 1][j] - 4 * T[i][j]);
            }
        }
    }

    // Function to compute R(u) on internal nodes of stagg-x mesh
    void ComputeRu(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Ru, double rho,
                   double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                Ru[i][j] = -rho * dx * (1.0 / 2 * (v[i][j - 1] + v[i][j]) * 1.0 / 2 * (u[i][j] + u[i - 1][j])
                                        - 1.0 / 2 * (v[i + 1][j - 1] + v[i + 1][j]) * 1.0 / 2 * (u[i][j] + u[i + 1][j])
                                        + 1.0 / 2 * (u[i][j] + u[i][j + 1]) * 1.0 / 2 * (u[i][j] + u[i][j + 1])
                                        - 1.0 / 2 * (u[i][j] + u[i][j - 1]) * 1.0 / 2 * (u[i][j] + u[i][j - 1]))
                           + mu * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]);
            }
        }
    }

    // Function to compute R(v) on internal nodes of stagg-y mesh
    void ComputeRv(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Rv, double rho,
                   double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rv[i][j] = -rho * dx * (1.0 / 2 * (v[i][j] + v[i - 1][j]) * 1.0 / 2 * (v[i][j] + v[i - 1][j])
                                        - 1.0 / 2 * (v[i][j] + v[i + 1][j]) * 1.0 / 2 * (v[i][j] + v[i + 1][j])
                                        + 1.0 / 2 * (u[i - 1][j + 1] + u[i][j + 1]) * 1.0 / 2 * (v[i][j] + v[i][j + 1])
                                        - 1.0 / 2 * (u[i - 1][j] + u[i][j]) * 1.0 / 2 * (v[i][j] + v[i][j - 1]))
                           + mu * (v[i - 1][j] + v[i + 1][j] + v[i][j - 1] + v[i][j + 1] - 4 * v[i][j]);
            }
        }
    }

    void Nusselt(int Nx, int Ny, double Tcold, double Thot, double L, double alpha, double dx, double &Nux_avg,
                 vector<vector<double> > &T1star, vector<vector<double> > &T1, vector<vector<double> > &qx,
                 vector<vector<double> > &u1star, vector<vector<double> > &u1, vector<vector<double> > &dTdx,
                 vector<double> &Nux, double dxstar) override {

        // Normalize distance, temperature and horizontal velocity
        dxstar = dx / L;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T1star[i][j] = (T1[i][j] - Tcold) / (Thot - Tcold);
            }
        }
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx + 1; j++) {
                u1star[i][j] = u1[i][j] * L / alpha;
            }
        }
        // Hot wall (x = 0) — dT/dx
        for (int i = 0; i < Ny; ++i) {
            dTdx[i][0] = (T1star[i][1] - T1star[i][0]) / dxstar;
            double T_central = 0.5 * (T1star[i][1] + T1star[i][0]);
            qx[i][0] = u1star[i][0] * T_central - dTdx[i][0];
        }
        // Internal cells (CDS)
        for (int i = 0; i < Ny; ++i) {
            for (int j = 1; j < Nx; ++j) {
                dTdx[i][j] = (T1star[i][j] - T1star[i][j - 1]) / dxstar;
                double T_central = 0.5 * (T1star[i][j] + T1star[i][j - 1]);
                qx[i][j] = u1star[i][j] * T_central - dTdx[i][j];
            }
        }
        // Vertical integration
        for (int j = 0; j < Nx; ++j) {
            double sum = 0.0;
            for (int i = 0; i < Ny; ++i) {
                sum += qx[i][j] * dxstar;
            }
            Nux[j] = sum;
        }
        // Average Nusselt
        Nux_avg = 0.0;
        for (int j = 0; j < Nx; ++j) {
            Nux_avg += Nux[j];
        }
        Nux_avg /= Nx;
    }
};

// Upwind Differencing Scheme (UDS) method
class UDSMethod : public ConvectiveScheme {

public:
    // Function to compute R(t) on each node of stagg-T mesh (same as stagg-P)
    void ComputeRt(double rho, int Ny, int Nx, double dx, double lambda, double cp, vector<vector<double> > &u,
                   vector<vector<double> > &v, vector<vector<double> > &Rt, vector<vector<double> > &T) override {
        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                double convW = u[i][j] > 0 ? u[i][j] * T[i][j - 1] : u[i][j] * T[i][j];

                double convE = u[i][j + 1] > 0 ? u[i][j + 1] * T[i][j] : u[i][j + 1] * T[i][j + 1];

                double convS = v[i + 1][j] > 0 ? v[i + 1][j] * T[i + 1][j] : v[i + 1][j] * T[i][j];

                double convN = v[i][j] > 0 ? v[i][j] * T[i][j] : v[i][j] * T[i - 1][j];

                Rt[i][j] = -rho * dx * cp * (convE - convW + convN - convS) +
                           lambda * (T[i][j + 1] + T[i][j - 1] + T[i - 1][j] + T[i + 1][j] - 4 * T[i][j]);
            }
        }
    }

    // Function to compute R(u) on internal nodes of stagg-x mesh
    void ComputeRu(vector<vector<double> > &u, vector<vector<double> > &v,
                   vector<vector<double> > &Ru, double rho, double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                double v_n = (v[i][j - 1] + v[i][j]) / 2;
                double conv_n = v_n > 0 ? v_n * u[i][j] : v_n * u[i - 1][j];

                double v_s = (v[i + 1][j - 1] + v[i + 1][j]) / 2;
                double conv_s = v_s > 0 ? v_s * u[i + 1][j] : v_s * u[i][j];

                double conv_w = u[i][j] > 0 ? u[i][j] * u[i][j - 1] : u[i][j] * u[i][j];

                double conv_e = u[i][j + 1] > 0 ? u[i][j + 1] * u[i][j] : u[i][j + 1] * u[i][j + 1];

                Ru[i][j] = -rho * dx * (conv_e - conv_w + conv_n - conv_s)
                           + mu * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]);
            }
        }
    }

    //Function to compute R(v) on internal nodes of stagg-y mesh
    void ComputeRv(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Rv, double rho,
                   double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                double conv_n = v[i][j] > 0 ? v[i][j] * v[i - 1][j] : v[i][j] * v[i][j];

                double conv_s = v[i + 1][j] > 0 ? v[i + 1][j] * v[i][j] : v[i + 1][j] * v[i + 1][j];

                double u_e = 0.5 * (u[i - 1][j + 1] + u[i][j + 1]);
                double conv_e = u_e > 0 ? u_e * v[i][j] : u_e * v[i][j + 1];

                double u_w = 0.5 * (u[i - 1][j] + u[i][j]);
                double conv_w = u_w > 0 ? u_w * v[i][j - 1] : u_w * v[i][j];

                Rv[i][j] = -rho * dx * (conv_e - conv_w + conv_n - conv_s)
                           + mu * (v[i - 1][j] + v[i + 1][j] + v[i][j - 1] + v[i][j + 1] - 4 * v[i][j]);
            }
        }
    }

    // Function to evaluate average Nusselt number
    void Nusselt(int Nx, int Ny, double Tcold, double Thot, double L, double alpha, double dx, double &Nux_avg,
                 vector<vector<double> > &T1star, vector<vector<double> > &T1, vector<vector<double> > &qx,
                 vector<vector<double> > &u1star, vector<vector<double> > &u1, vector<vector<double> > &dTdx,
                 vector<double> &Nux, double dxstar) override {
        // Normalize distance, temperature and horizontal velocity
        dxstar = dx / L;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T1star[i][j] = (T1[i][j] - Tcold) / (Thot - Tcold);
            }
        }
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx + 1; j++) {
                u1star[i][j] = u1[i][j] * L / alpha;
            }
        }
        // Heat flux (UDS)
        for (int i = 0; i < Ny; ++i) {
            // (x = 0)
            dTdx[i][0] = (T1star[i][1] - T1star[i][0]) / dxstar;
            double T_upwind_0 = (u1star[i][0] > 0.0) ? T1star[i][0] : T1star[i][1];
            qx[i][0] = u1star[i][0] * T_upwind_0 - dTdx[i][0];

            // (j = 1)
            dTdx[i][1] = (T1star[i][1] - T1star[i][0]) / dxstar;
            double T_upwind_1 = (u1star[i][1] > 0.0) ? T1star[i][0] : T1star[i][1];
            qx[i][1] = u1star[i][1] * T_upwind_1 - dTdx[i][1];

            // Internals (j ≥ 2)
            for (int j = 2; j < Nx; ++j) {
                if (u1star[i][j] > 0.0) {
                    dTdx[i][j] = (T1star[i][j - 1] - T1star[i][j - 2]) / dxstar;
                } else {
                    dTdx[i][j] = (T1star[i][j] - T1star[i][j - 1]) / dxstar;
                }
                double T_upwind = (u1star[i][j] > 0.0) ? T1star[i][j - 1] : T1star[i][j];
                qx[i][j] = u1star[i][j] * T_upwind - dTdx[i][j];
            }
        }
        // Vertical integration
        for (int j = 0; j < Nx; ++j) {
            double sum = 0.0;
            for (int i = 0; i < Ny; ++i) {
                sum += qx[i][j] * dxstar;
            }
            Nux[j] = sum;
        }
        // Average Nusselt
        Nux_avg = 0.0;
        for (int j = 0; j < Nx; ++j) {
            Nux_avg += Nux[j];
        }
        Nux_avg /= Nx;
    }
};

// Function to choose the method
ConvectiveScheme* createCalculator(const string &method) {
    if (method == "CDS") {
        return new CDSMethod();
    } else if (method == "UDS") {
        return new UDSMethod();
    } else {
        cerr << "Not a valid method!" << endl;
        return nullptr;
    }
}

// Function to apply boundary conditions
void BoundaryConditions(int Nx, int Ny, double Thot, double Tcold, vector<vector<double> > &Tn_1,
                        vector<vector<double> > &Tn, vector<vector<double> > &T1) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            // LEFT WALL (x = 0): Dirichlet BC
            if (j == 0) {
                Tn_1[i][j] = Thot;
                Tn[i][j] = Thot;
                T1[i][j] = Thot;
            }
            // RIGHT WALL (x = L): Dirichlet BC
            else if (j == Nx - 1) {
                Tn_1[i][j] = Tcold;
                Tn[i][j] = Tcold;
                T1[i][j] = Tcold;
            }
        }
    }
    // TOP WALL (y = L -> i = 0): Neumann dT/dz = 0 -> T[0][j] = T[1][j]
    for (int j = 0; j < Nx; j++) {
        Tn_1[0][j] = Tn_1[1][j];
        Tn[0][j] = Tn[1][j];
        T1[0][j] = T1[1][j];
    }
    // BOTTOM WALL (y = 0 -> i = Ny - 1): Neumann dT/dz = 0 -> T[Ny-1][j] = T[Ny-2][j]
    for (int j = 0; j < Nx; j++) {
        Tn_1[Ny - 1][j] = Tn_1[Ny - 2][j];
        Tn[Ny - 1][j] = Tn[Ny - 2][j];
        T1[Ny - 1][j] = T1[Ny - 2][j];
    }
}

// Function to solve the pressure field with Gauss-Seidel solver on all nodes of stagg-P mesh
void PoissonSolver(double maxResP, double maxIteP, double rho, double dx, double dt, int Nx, int Ny, double omega_P,
                   vector<vector<double> > &P1, vector<vector<double> > &vP, vector<vector<double> > &uP,
                   vector<vector<double> > &Pg) {
    double resP = maxResP + 1;
    int iteP = 0;
    while (resP > maxResP && iteP < maxIteP) {
        double maxDiffP = 0.0;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (i == 0)
                    P1[i][j] = 1.0 / 1 * (P1[i+1][j] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (j == 0)
                    P1[i][j] = 1.0 / 1 * (P1[i][j+1] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (i == Ny - 1)
                    P1[i][j] = 1.0 / 1 * (P1[i-1][j] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (j == Nx - 1)
                    P1[i][j] = 1.0 / 1 * (P1[i][j-1] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else
                    P1[i][j] = 1.0 / 4 * (P1[i+1][j] + P1[i][j+1] + P1[i-1][j] + P1[i][j-1]
                        - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                double diffP = fabs(P1[i][j] - Pg[i][j]);
                if (diffP > maxDiffP) {
                    maxDiffP = diffP;
                }
            }
        }
        resP = maxDiffP;
        // Update P guess with over-relaxation
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (omega_P < 1e-8) {
                    Pg[i][j] = P1[i][j];
                } else {
                    Pg[i][j] = Pg[i][j] + omega_P * (P1[i][j] - Pg[i][j]);
                }
            }
        }
        iteP++;
    }
}

int main () {
    // Start timer
    auto start = high_resolution_clock::now();

    // Physical data
    int N;
    cout << "Insert number of control volumes: ";                                   // Ask for mesh refinement
    cin >> N;
    int Nx = N;
    int Ny = N;
    double Ra;
    double L = 1.0;
    double dx = L / Nx;
    double Pr = 0.71;
    double Thot = 1.0;
    double Tcold = 0.0;
    double Tinf = 0.5;
    double cp = 0.71;
    double lambda = 1.0;
    double mu = Pr * lambda / cp;
    double rho = 1.0;
    double beta = 1.0;
    cout << "Insert Rayleigh number: ";                                       // Ask for desired Rayleigh number
    cin >> Ra;
    double g = Ra * mu * lambda / (cp * beta * pow(rho, 2) * (Thot - Tcold) * pow(L, 3));
    double alpha = lambda / (rho * cp);
    double Nux_avg = 0.0;
    double Vmax;
    double dxstar;

    // Numerical data
    double maxResP = 1e-6;
    double maxIte = 1e6;
    double maxResT = 1e-4;
    double maxResU = maxResT;
    double maxResV = maxResT;
    double t_count = 0.0;
    double res1 = maxResU + 1;
    double res2 = maxResV + 1;
    double res3 = maxResT + 1;
    double dt = 1e-4;
    double dtc, dtd, dtt;
    int ite = 0;
    string method;
    string graph;
    double omega_T = 0.5;                                                     // Relaxation factors
    double omega_u = 0.5;
    double omega_v = 0.5;
    double omega_P = 1.3;

    cout << "Choose the method ('CDS', 'UDS', or 'AUTO'): ";                 // Ask for convective scheme
    cin >> method;
    ranges::transform(method, method.begin(), ::toupper);

    ConvectiveScheme* calculator = nullptr;
    ConvectiveScheme* calculator_UDS = nullptr;
    ConvectiveScheme* calculator_CDS = nullptr;

    bool switched = false;
    double res3_old = res3;

    if (method == "CDS" || method == "UDS") {
        calculator = createCalculator(method);
        if (!calculator) {
            cerr << "Error: Invalid method selected." << endl;
            return 1;
        }
        cout << "Running with fixed " << method << " scheme.\n";
    } else if (method == "AUTO") {
        calculator_UDS = createCalculator("UDS");
        calculator_CDS = createCalculator("CDS");
        if (!calculator_UDS || !calculator_CDS) {
            cerr << "Error: Could not create convection schemes." << endl;
            return 1;
        }
        calculator = calculator_UDS; // start with UDS
        cout << "AUTO mode: starting with UDS, will switch to CDS based on residuals.\n";
    } else {
        cerr << "Error: Method must be 'CDS', 'UDS', or 'AUTO'." << endl;
        return 1;
    }

    // Pressure, Temperature, u and v matrices
    auto P1    = CreateMatrix(Ny, Nx);
    auto Pg    = CreateMatrix(Ny, Nx);
    auto Tn_1 = CreateMatrix(Ny, Nx);
    auto Tn = CreateMatrix(Ny, Nx);
    auto T1 = CreateMatrix(Ny, Nx);
    auto un_1  = CreateMatrix(Ny, Nx + 1);
    auto un    = CreateMatrix(Ny, Nx + 1);
    auto u1    = CreateMatrix(Ny, Nx + 1);
    auto uP    = CreateMatrix(Ny, Nx + 1);
    auto vn_1  = CreateMatrix(Ny + 1, Nx);
    auto vn    = CreateMatrix(Ny + 1, Nx);
    auto v1    = CreateMatrix(Ny + 1, Nx);
    auto vP    = CreateMatrix(Ny + 1, Nx);

    // Convective-diffusive, Buoyancy term R() matrices
    auto Ru1    = CreateMatrix(Ny, Nx + 1);
    auto Run    = CreateMatrix(Ny, Nx + 1);
    auto Run_1  = CreateMatrix(Ny, Nx + 1);
    auto Rv1    = CreateMatrix(Ny + 1, Nx);
    auto Rvn    = CreateMatrix(Ny + 1, Nx);
    auto Rvn_1  = CreateMatrix(Ny + 1, Nx);
    auto Rt1 = CreateMatrix(Ny, Nx);
    auto Rtn = CreateMatrix(Ny, Nx);
    auto Rtn_1 = CreateMatrix(Ny, Nx);
    auto Rbn = CreateMatrix(Ny + 1, Nx);
    auto Rbn_1 = CreateMatrix(Ny + 1, Nx);

    // Adimensional analysis matrices and vectors
    auto T1star = CreateMatrix(Ny, Nx);
    auto u1star = CreateMatrix(Ny, Nx + 1);
    auto v1star = CreateMatrix(Ny + 1, Nx);
    auto dTdx = CreateMatrix(Ny, Nx);
    auto qx = CreateMatrix(Ny, Nx);
    auto V = CreateMatrix(Ny, Nx);
    vector Nux(Ny, 0.0);
    vector u1star_L2(Ny, 0.0);
    vector v1star_L2(Ny, 0.0);

    // Initial velocity fields = 0.0
    for (int i = 1; i < Ny; i++) {
        for (int j = 1; j < Nx + 1; j++) {
            u1[i][j] = 0.0;
            un[i][j] = u1[i][j];
            un_1[i][j] = un[i][j];
        }
    }
    for (int i = 1; i < Ny + 1; i++) {
        for (int j = 1; j < Nx; j++) {
            v1[i][j] = 0.0;
            vn[i][j] = v1[i][j];
            vn_1[i][j] = vn[i][j];
        }
    }
    // Initial pressure field = 0.0 and linearly distributed temperature field
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            Pg[i][j] = 0.0;
            T1[i][j] = 1.0 - i * dx;
            Tn[i][j] = T1[i][j];
            Tn_1[i][j] = Tn[i][j];
        }
    }

    // Apply temperature boundary conditions
    BoundaryConditions(Nx, Ny, Thot, Tcold, Tn_1, Tn, T1);

    // Compute R()^n-1
    calculator->ComputeRu(un_1, vn_1, Run_1, rho, dx, mu, Nx, Ny);
    calculator->ComputeRv(un_1, vn_1, Rvn_1, rho, dx, mu, Nx, Ny);
    calculator->ComputeRt(rho, Ny, Nx, dx, lambda, cp, un_1, vn_1, Rtn_1, Tn_1);

    // Compute R()^n
    calculator->ComputeRu(un, vn, Run, rho, dx, mu, Nx, Ny);
    calculator->ComputeRv(un, vn, Rvn, rho, dx, mu, Nx, Ny);
    calculator->ComputeRt(rho, Ny, Nx, dx, lambda, cp, un, vn, Rtn, Tn);

    // Time loop until convergence is reached
    while ((res1 > maxResU || res2 > maxResV || res3 > maxResT) && ite < maxIte) {
        double maxDiff1 = 0.0;
        double maxDiff2 = 0.0;
        double maxDiff3 = 0.0;

        // Step 1a
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                if (ite < 5) {                                         // Using Euler Explicit for the first time steps
                    uP[i][j] = un[i][j] + dt / (rho * pow(dx, 2)) * Run[i][j];
                  } else {
                      uP[i][j] = un[i][j] + dt / (rho * pow(dx, 2)) * (3.0 / 2 * Run[i][j] - 1.0 / 2 * Run_1[i][j]);
                  }
            }
        }

        // Step 1b + Boussinesq approx.
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rbn[i][j] = rho * pow(dx, 2) * beta * (Tn[i][j] - Tinf) * g;
                Rbn_1[i][j] = rho * pow(dx, 2) * beta * (Tn_1[i][j] - Tinf) * g;
                if (ite < 5) {
                    vP[i][j] = vn[i][j] + dt / (rho * pow(dx, 2)) * (Rvn[i][j] + Rbn[i][j]);
                } else {
                    vP[i][j] = vn[i][j] + dt / (rho * pow(dx, 2)) * (
                                   3.0 / 2 * Rvn[i][j] - 1.0 / 2 * Rvn_1[i][j] + 3.0 / 2 * Rbn[i][j] - 1.0 / 2 * Rbn_1[
                                       i][j]);
                }
            }
        }

        // Step 2
        PoissonSolver(maxResP, maxIte, rho, dx, dt, Nx, Ny, omega_P, P1, vP, uP, Pg);

        // Step 3
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                u1[i][j] = uP[i][j] - dt / (rho * dx) * (P1[i][j] - P1[i][j - 1]);
                double diff1 = fabs(u1[i][j] - un[i][j]);
                if (diff1 > maxDiff1) {
                    maxDiff1 = diff1;
                }
            }
        }
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                v1[i][j] = vP[i][j] - dt / (rho * dx) * (P1[i - 1][j] - P1[i][j]);
                double diff2 = fabs(v1[i][j] - vn[i][j]);
                if (diff2 > maxDiff2) {
                    maxDiff2 = diff2;
                }
            }
        }

        // Temperature boundary conditions
        BoundaryConditions(Nx, Ny, Thot, Tcold, Tn_1, Tn, T1);

        // Temperature evaluation
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                if (ite < 5) {
                    T1[i][j] = Tn[i][j] + dt / (pow(dx, 2) * rho * cp) * Rtn[i][j];
                } else {
                    T1[i][j] = Tn[i][j] + dt / (pow(dx, 2) * rho * cp) * (3.0 / 2 * Rtn[i][j] - 1.0 / 2 * Rtn_1[i][j]);
                }
                double diff3 = fabs(T1[i][j] - Tn[i][j]);
                if (diff3 > maxDiff3) {
                    maxDiff3 = diff3;
                }
            }
        }

        // Update residuals
        t_count += dt;
        res1 = maxDiff1;
        res2 = maxDiff2;
        res3 = maxDiff3;
        // Exits the loop if res3 < maxResT; if not, continue

        if (method == "AUTO" && !switched && res3 < 1e-2 && fabs(res3 - res3_old) / res3_old < 5e-2) {
            calculator = calculator_CDS;
            switched = true;
            cout << "Switched from UDS to CDS at iteration " << ite << " (res3 = " << res3 << ")\n";
        }
        res3_old = res3;

        // Under-relaxation and update u^n and u^n-1
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                if (omega_u < 1e-8) {
                    un_1[i][j] = un[i][j];
                    un[i][j] = u1[i][j];
                } else {
                    u1[i][j] = un[i][j] + omega_u * (u1[i][j] - un[i][j]);
                    un_1[i][j] = un[i][j];
                    un[i][j] = u1[i][j];
                }
            }
        }

        // Under-relaxation and update v^n and v^n-1
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                if (omega_v < 1e-8) {
                    vn_1[i][j] = vn[i][j];
                    vn[i][j] = v1[i][j];
                } else {
                    v1[i][j] = vn[i][j] + omega_v * (v1[i][j] - vn[i][j]);
                    vn_1[i][j] = vn[i][j];
                    vn[i][j] = v1[i][j];
                }
            }
        }

        // Under-relaxation and update T^n and T^n-1
        if (omega_T < 1e-8) {
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    Tn_1[i][j] = Tn[i][j];
                    Tn[i][j] = T1[i][j];
                }
            }
        } else {
            for (int i = 1; i < Ny - 1; i++) {
                for (int j = 1; j < Nx - 1; j++) {
                    T1[i][j] = Tn[i][j] + omega_T * (T1[i][j] - Tn[i][j]);
                }
            }
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    Tn_1[i][j] = Tn[i][j];
                    Tn[i][j] = T1[i][j];
                }
            }
        }

        // Compute R()^n+1
        calculator->ComputeRu(u1, v1, Ru1, rho, dx, mu, Nx, Ny);
        calculator->ComputeRv(u1, v1, Rv1, rho, dx, mu, Nx, Ny);
        calculator->ComputeRt(rho, Ny, Nx, dx, lambda, cp, u1, v1, Rt1, T1);

        // Update R()^n and R()^n-1
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                Run_1[i][j] = Run[i][j];
                Run[i][j] = Ru1[i][j];
            }
        }
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rvn_1[i][j] = Rvn[i][j];
                Rvn[i][j] = Rv1[i][j];
            }
        }
        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                Rtn_1[i][j] = Rtn[i][j];
                Rtn[i][j] = Rt1[i][j];
            }
        }

        //Update time step
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                V[i][j] = sqrt(pow(u1[i][j], 2) + pow(v1[i][j], 2));
            }
        }
        Vmax = std::ranges::max(std::views::join(V));
        dtc = 0.35 * dx * L / Vmax;
        dtd = 0.20 * pow(dx, 2) / (mu / rho);
        dtt = 0.20 * pow(dx, 2) / (lambda / (rho * cp));
        dt = min ({dtc, dtd, dtt});

        // Go to next time step
        ite++;
        cout << "ite = " << ite << " | res1 = " << res1 << " | res2 = " << res2 << " | res3 = " << res3 << " | dt = " << dt << endl;
    }

    cout << "Steady state reached in " << t_count << " seconds" << endl;

    // .txt files for plotting

    // File to print temperature distribution for Plot 1
    ofstream TemperatureDistribution("TemperatureDistribution.txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            TemperatureDistribution << j * dx + dx / 2 << " " << L - i * dx - dx / 2 << " " << T1[i][j] << endl;
        }
        TemperatureDistribution << "\n";
    }
    TemperatureDistribution.close();

    // File to print velocity u distribution for Plot 2
    ofstream VelocityUDistribution("VelocityUDistribution.txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx + 1; j++) {
            VelocityUDistribution << j * dx << " " << L - i * dx - dx / 2 << " " << u1[i][j] << endl;
        }
        VelocityUDistribution << "\n";
    }
    VelocityUDistribution.close();

    // File to print velocity v distribution for Plot 3
    ofstream VelocityVDistribution("VelocityVDistribution.txt");
    for (int i = 0; i < Ny + 1; i++) {
        for (int j = 0; j < Nx; j++) {
            VelocityVDistribution << j * dx + dx / 2 << " " << L - i * dx << " " << v1[i][j] << endl;
        }
        VelocityVDistribution << "\n";
    }
    VelocityVDistribution.close();

    // File to print u(y) at L/2 for Plot 4
    ofstream Uyl2("Uy.txt");
    Uyl2 << 1 << " " << 0 << endl;
    for (int i = 0; i < Ny; i++) {
        Uyl2 << L - i * dx - dx / 2 << " " << (u1[i][Nx / 2] + u1[i][Nx / 2 + 1]) / 2 << endl;
    }
    Uyl2 << 0 << " " << 0 << endl;
    Uyl2.close();

    // File to print v(x) at L/2 for Plot 5
    ofstream Vxl2("Vx.txt");
    Vxl2 << 0 << " " << 0 << endl;
    for (int j = 0; j < Nx; j++) {
        Vxl2 << j * dx + dx / 2 << " " << (v1[Ny / 2][j] + v1[Ny / 2 + 1][j]) / 2 << endl;
    }
    Vxl2 << 1 << " " << 0 << endl;
    Vxl2.close();

    // Average Nusselt number
    calculator->Nusselt(Nx, Ny, Tcold, Thot, L, alpha, dx, Nux_avg, T1star, T1, qx, u1star, u1, dTdx, Nux, dxstar);
    cout << "Average Nusselt for Ra = " << Ra << " --> " << Nux_avg << endl;

    //Maximum u* velocity at x = L/2
    for (int i = 0; i < Ny; i++) {
        u1star_L2[i] = (u1[i][Nx / 2] + u1[i][Nx / 2 + 1]) / 2 * L / alpha;
    }
    auto ustar_max = ranges::max_element(u1star_L2);
    int u_index = distance(u1star_L2.begin(), ustar_max);
    double y_u = L - (u_index + 0.5) * dx;
    double ystar_u = y_u / L;
    cout << "Max u* velocity = " << *ustar_max << " at y/L = " << ystar_u << endl;

    // Maximum v* velocity at y = L/2
    for (int j = 0; j < Nx; j++) {
        v1star_L2[j] = (v1[Ny / 2][j] + v1[Ny / 2 + 1][j]) / 2 * L / alpha;
    }
    auto vstar_max = ranges::max_element(v1star_L2);
    int v_index = distance(v1star_L2.begin(), vstar_max);
    double x_v = (v_index + 0.5) * dx;
    double xstar_v = x_v / L;
    cout << "Max v* velocity = " << *vstar_max << " at x/L = " << xstar_v << endl;

    // Stop timer and print total duration
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    cout << "Total Execution Time: " << duration.count() << " seconds" << endl;

    if (method == "AUTO") {
        delete calculator_UDS;
        delete calculator_CDS;
    } else {
        delete calculator;
    }

    return 0;
}