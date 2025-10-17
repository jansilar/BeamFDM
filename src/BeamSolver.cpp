#include "BeamSolver.h"
#include <Eigen/Dense>

BeamSolver::BeamSolver(int N_, double L_, double E_, double I_, double rho_, double A_, double c_)
    : N(N_), L(L_), E(E_), I(I_), rho(rho_), A(A_), c(c_)
{
    dx = L / (N - 1);
    x.resize(N);
    y.resize(N, 0.0);
    yPrev.resize(N, 0.0);
    yNew.resize(N, 0.0);
    for (int i = 0; i < N; ++i)
        x[i] = i * dx;
}

void BeamSolver::solveStatic(double q) {
    // q .. rozložené zatížení [N/m]
    // Bernoulli-Navier eq: EI * D4(y) = q
    // pětidiagonální systém: D4(y) = (y[i-2] - 4y[i-1] + 6y[i] - 4y[i+1] + y[i+2]) / h^4
    // pro i = 2..N-3 (vnitřní uzly)
    // BC: y(0)=0, y'(0)=0, y''(L)=0, y'''(L)=0 (clamped-free beam)
    // systém: A*y = f

    using namespace Eigen;
    MatrixXd A = MatrixXd::Zero(N, N);
    VectorXd f = VectorXd::Zero(N);

    // vnitřní uzly (4. derivace centrálně)
    for (int i = 2; i < N - 2; ++i) {
        A(i,i - 2) =  1.0;
        A(i,i - 1) = -4.0;
        A(i,i)     =  6.0;
        A(i,i + 1) = -4.0;
        A(i,i + 2) =  1.0;
        f(i) = q*std::pow(dx, 4) / (E * I);
    }

    // ---- Okrajové podmínky (clamped-free beam) ----
    // y(0)=0, y'(0)=0, y''(L)=0, y'''(L)=0
    // Implementace: nahradíme první a poslední 2 řádky rovnicemi pro BC

    // y(0)=0
    A(0,0) = 1.0; f(0) = 0.0;
    // y'(0)=0 → (-3y0 + 4y1 - y2) / (2h) = 0
    A(1,0) = -3.0;
    A(1,1) =  4.0;
    A(1,2) = -1.0;
    f(1) = 0.0;

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

    // ---- Vyřeš lineární systém ----
    
    Eigen::VectorXd y_vec = A.fullPivLu().solve(f);
    for (int i = 0; i < N; ++i)
        y[i] = y_vec(i);
}

void BeamSolver::stepDynamic(double q, double dt, int steps) {
    //q .. rozložené zatížení [N/m]
    //dt .. časový krok [s]
    //steps .. počet časových kroků
    double M = rho*A/(dt*dt);       // hmotnostní člen
    double D = c/(2*dt);            // tlumící člen
    double K = E*I/std::pow(dx,4);  // tuhostní člen
    yNew[0] = 0.0; // BC y(0)=0
    yNew[1] = 0.0; // BC y'(0)=0
    for (int step = 0; step < steps; ++step) {
        // Vnitřní uzly
        for (int i = 2; i < N - 2; ++i) {
            double D4y = (y[i-2] - 4*y[i-1] + 6*y[i] - 4*y[i+1] + y[i+2]);
            double f = q; // externí síla
            yNew[i] = (f + M*(2*y[i] - yPrev[i]) + D*yPrev[i] - K*D4y) / (M + D);
        }
        // BC na volném konci
        // y''(L)=0 → (yN-3 - 2yN-2 + yN-1)/h^2 = 0
        // y'''(L)=0 → (-yN-4 + 4yN-3 -5yN-2 + 2yN-1)/h^3 = 0
        // vyřešit pro neznámé yN-1, yN-2
        // z první: yN-1 = 2yN-2 - yN-3
        // dosadit do druhé: -yN-4 + 4yN-3 -5yN-2 + 2(2yN-2 - yN-3) = 0
        // -yN-4 + 2yN-3 - yN-2 = 0
        // yN-2 = 2yN-3 - yN-4
        // dosadit zpět: yN-1 = 2(2yN-3 - yN-4) - yN-3 = 3yN-3 - 2yN-4
        // tedy:
        yNew[N-2] = 2*yNew[N-3] - yNew[N-4];
        yNew[N-1] = 3*yNew[N-3] - 2*yNew[N-4];

        // Posun v čase
        yPrev = y;
        y = yNew;
    }
}

double BeamSolver::analyticMaxDeflection(double q) const {
    // analytické maximum pro konzolu s uniform q: y_max = qL^4 / (8EI)
    return q * std::pow(L, 4) / (8.0 * E * I);
}