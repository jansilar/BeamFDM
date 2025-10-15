#include "BeamSolver.h"

BeamSolver::BeamSolver(int N_, double L_, double E_, double I_, double q_)
    : N(N_), L(L_), E(E_), I(I_), q(q_)
{
    h = L / (N - 1);
    x.resize(N);
    y.resize(N, 0.0);
    for (int i = 0; i < N; ++i)
        x[i] = i * h;
}

void BeamSolver::solve() {
    // pětidiagonální systém: EI * D4(y) = q
    // D4(y) = (y[i-2] - 4y[i-1] + 6y[i] - 4y[i+1] + y[i+2]) / h^4
    // pro i = 2..N-3 (vnitřní uzly)
    // BC: y(0)=0, y'(0)=0, y''(L)=0, y'''(L)=0 (clamped-free beam)
    // využijeme metodu Gaussovy eliminace pro pásovou matici
    // (zde malý N → klasická implementace OK)

    // systém: A*y = f
    std::vector<std::vector<double>> A(N, std::vector<double>(N, 0.0));
    std::vector<double> f(N, 0); // uniform q

    // vnitřní uzly (4. derivace centrálně)
    for (int i = 2; i < N - 2; ++i) {
        A[i][i - 2] =  1.0;
        A[i][i - 1] = -4.0;
        A[i][i]     =  6.0;
        A[i][i + 1] = -4.0;
        A[i][i + 2] =  1.0;
        f[i] = q*std::pow(h, 4) / (E * I);
    }

    // ---- Okrajové podmínky (clamped-free beam) ----
    // y(0)=0, y'(0)=0, y''(L)=0, y'''(L)=0
    // Implementace: nahradíme první a poslední 2 řádky rovnicemi pro BC

    // y(0)=0
    A[0][0] = 1.0; f[0] = 0.0;
    // y'(0)=0 → (-3y0 + 4y1 - y2) / (2h) = 0
    A[1][0] = -3.0/(2*h);
    A[1][1] =  4.0/(2*h);
    A[1][2] = -1.0/(2*h);
    f[1] = 0.0;

    // y''(L)=0 → (yN-3 - 2yN-2 + yN-1)/h^2 = 0
    A[N-2][N-3] =  1.0/(h*h);
    A[N-2][N-2] = -2.0/(h*h);
    A[N-2][N-1] =  1.0/(h*h);
    f[N-2] = 0.0;

    // y'''(L)=0 → (-yN-4 + 4yN-3 -5yN-2 + 2yN-1)/h^3 = 0
    A[N-1][N-4] = -1.0/(h*h*h);
    A[N-1][N-3] =  4.0/(h*h*h);
    A[N-1][N-2] = -5.0/(h*h*h);
    A[N-1][N-1] =  2.0/(h*h*h);
    f[N-1] = 0.0;

    // ---- Vyřeš lineární systém ----
    // jednoduchá Gaussova eliminace
    for (int k = 0; k < N; ++k) {
        // pivot
        double piv = A[k][k];
        if (std::fabs(piv) < 1e-14) continue;
        for (int j = k; j < N; ++j)
            A[k][j] /= piv;
        f[k] /= piv;

        // eliminace
        for (int i = k + 1; i < N; ++i) {
            double factor = A[i][k];
            for (int j = k; j < N; ++j)
                A[i][j] -= factor * A[k][j];
            f[i] -= factor * f[k];
        }
    }

    // zpětný průchod
    for (int i = N - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < N; ++j)
            sum += A[i][j] * y[j];
        y[i] = f[i] - sum;
    }
}

double BeamSolver::analyticMaxDeflection() const {
    // analytické maximum pro konzolu s uniform q: y_max = qL^4 / (8EI)
    return q * std::pow(L, 4) / (8.0 * E * I);
}