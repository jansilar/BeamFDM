#include "BeamSolver.h"
#include <Eigen/Dense>
#include <algorithm>

BeamSolver::BeamSolver(int N_, double L_, double E_, double I_, double rho_, double A_, double c_)
    : N(N_), L(L_), E(E_), I(I_), rho(rho_), A(A_), c(c_)
{
    if (N < 5)
        throw std::invalid_argument("N must be at least 5 for the beam static solver.");
    dx = L / (N - 1);
    x.resize(N);
    y.resize(N, 0.0);
    yPrev.resize(N, 0.0);
    yNew.resize(N, 0.0);
    for (int i = 0; i < N; ++i)
        x[i] = i * dx;
}

void BeamSolver::solveStatic(double q) {
    // q .. distrubuted load [N/m]
    // Bernoulli-Navier eq: EI * D4(y) = q
    // pentadiagonal system (stencil): D4(y) = (y[i-2] - 4y[i-1] + 6y[i] - 4y[i+1] + y[i+2]) / h^4
    // for i = 2..N-3 (internal nodes)
    // BC: y(0)=0, y'(0)=0, y''(L)=0, y'''(L)=0 (clamped-free beam)
    // system: A*y = f

    using namespace Eigen;
    MatrixXd A = MatrixXd::Zero(N, N);
    VectorXd f = VectorXd::Zero(N);

    // inner nodes - 4th derivaive - central finite difference
    for (int i = 2; i < N - 2; ++i) {
        A(i,i - 2) =  1.0;
        A(i,i - 1) = -4.0;
        A(i,i)     =  6.0;
        A(i,i + 1) = -4.0;
        A(i,i + 2) =  1.0;
        f(i) = q*std::pow(dx, 4) / (E * I);
    }

    // ---- BC (clamped-free beam) ----
    // y(0)=0, y'(0)=0, y''(L)=0, y'''(L)=0

    // First two rows of A

    // y(0)=0
    A(0,0) = 1.0; f(0) = 0.0;
    // y'(0)=0 → (-3y0 + 4y1 - y2) / (2h) = 0
    A(1,0) = -3.0;
    A(1,1) =  4.0;
    A(1,2) = -1.0;
    f(1) = 0.0;

    //Last two rows of A

    // y''(L)=0 → (yN-3 - 2yN-2 + yN-1)/h^2 = 0
    A(N-2,N-3) =  1.0;
    A(N-2,N-2) = -2.0;
    A(N-2,N-1) =  1.0;
    f(N-2) = 0.0;

    // y'''(L)=0 → (-yN-4 + 4yN-3 -5yN-2 + 2yN-1)/h^3 = 0
    A(N-1,N-4) = -1.0;
    A(N-1,N-3) =  4.0;
    A(N-1,N-2) = -5.0;
    A(N-1,N-1) =  2.0;
    f(N-1) = 0.0;

    // ---- Solve the linear system ----
    
    Eigen::VectorXd y_vec = A.fullPivLu().solve(f);
    for (int i = 0; i < N; ++i)
        y[i] = y_vec(i);
}

void BeamSolver::stepDynamic(double q, double dt, int steps) {
    //q .. distributed load [N/m]
    //dt .. time step [s]
    //steps .. number of time steps to perform
    // Solves the equation Bernoulli-Navier beam equation with
    // damping and inertia terms:
    // rho A d^2y/dt^2 + c dy/dt} + E I d^4y/dx^4} = q
    // using explicit time stepping, central differences in time and space
    // d^2y/dt^2 ≈ (y_new - 2y + y_prev) / dt^2
    // dy/dt ≈ (y_new - y_prev) / (2 dt)
    // d^4y/dx^4 ≈ (y[i-2] - 4y[i-1] + 6y[i] - 4y[i+1] + y[i+2]) / dx^4
    double M = rho*A/(dt*dt);       // mass therm
    double D = c/(2*dt);            // damping term
    double K = E*I/std::pow(dx,4);  // stiffness term
    // Clamped end BC:
    yNew[0] = 0.0; // BC y(0)=0
    yNew[1] = 0.0; // BC y'(0)=0
    for (int step = 0; step < steps; ++step) {
        // Inner nodes
        for (int i = 2; i < N - 2; ++i) {
            //fourth derivative stencil
            double D4y = (y[i-2] - 4*y[i-1] + 6*y[i] - 4*y[i+1] + y[i+2]);
            // Update state - equation derived from the PDE discretization
            yNew[i] = (q + M*(2*y[i] - yPrev[i]) + D*yPrev[i] - K*D4y) / (M + D);
        }
        // Free end BCs:
        // y''(L)=0 → (yN-3 - 2yN-2 + yN-1)/h^2 = 0
        // y'''(L)=0 → (-yN-4 + 4yN-3 -5yN-2 + 2yN-1)/h^3 = 0
        // solve for yN-1, yN-2
        yNew[N-2] = 2*yNew[N-3] - yNew[N-4];
        yNew[N-1] = 3*yNew[N-3] - 2*yNew[N-4];

        // Advance time
        std::swap(yPrev,y);
        std::swap(y, yNew);
    }
}

void BeamSolver::resetState() {
    // Reset state to initial conditions (all zeros)
    for (int i = 0; i < N; ++i) {
        y[i] = 0.0;
        yPrev[i] = 0.0;
        yNew[i] = 0.0;
    }
}   


double BeamSolver::analyticMaxDeflection(double q) const {
    // max deflection using analytical solution of static problem
    return q * std::pow(L, 4) / (8.0 * E * I);
}