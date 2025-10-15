#include "BeamSolver.h"
#include <iostream>
#include <iomanip>
#include <algorithm> 

int main() {
    // Parametry nosníku
    double L = 2.0;           // délka [m]
    double E = 1.0e10;        // Young [Pa]
    double d = 0.10;          // průměr [m]
    double I = M_PI * std::pow(d, 4) / 64.0; // moment setrvačnosti [m^4]
    double q = 300.0;         // zatížení [N/m]

    BeamSolver solver(30, L, E, I, q);
    solver.solve();

    const auto& y = solver.getY();
    const auto& x = solver.getX();

    std::cout << "x [m]\t y [mm]\n";
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << std::fixed << std::setprecision(3)
                  << x[i] << "\t" << y[i]*1000 << "\n";
    }

    double yAnal = solver.analyticMaxDeflection();
    std::cout << "\nAnalytic max deflection = " << yAnal*1000 << " mm\n";
    std::cout << "Numerical max deflection = " << (*std::max_element(y.begin(), y.end()))*1000 << " mm\n";
}