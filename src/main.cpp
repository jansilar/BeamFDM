#include "BeamSolver.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>

int main() {
    int N = 10;          // počet uzlů
    std::string fileName = "beam_results.csv";
    // Parametry prutu
    double L = 2.0;           // délka [m]
    double E = 1.0e10;        // Young [Pa]
    double d = 0.10;          // průměr [m]
    double I = M_PI * std::pow(d, 4) / 64.0; // moment setrvačnosti [m^4]
    double q = 300.0;         // zatížení [N/m]
    // Parametry pro dynamiku
    double rho = 700.0;      // hustota [kg/m^3]
    double A = M_PI * std::pow(d, 2) / 4.0; // průřezová plocha [m^2]
    double omega1 = std::pow(1.875,2) * std::sqrt(E*I/(rho*A*std::pow(L,4))); // 1. vlastní frekvence vetknutého nosníku
    double ksi = 0.05;        // poměr krytického tlumení (2%)
    double c = ksi*(2*omega1*rho*A);        // tlumení
    double dx = L / (N - 1);          // krok sítě [m]
    double dt = 0.01*std::pow(dx,2)/2*std::sqrt(rho*A/(E*I));        // časový krok [s]

    BeamSolver solver(N, L, E, I, q, rho, A, c);
    const auto& x = solver.getX();
    const auto& y = solver.getY();

    // Start profiling
    auto start = std::chrono::high_resolution_clock::now();
    // Run solver
    //static:
    //solver.solveStatic();
    //dynamic:
    double simTime = 1; // celkový čas simulace [s]
    double t = 0;
    while (t < simTime){
        solver.stepDynamic(dt, 1);
        t += dt;
        std::cout << "t = " << t << " s, yLast = " << y[N-1]*1000 << " mm\n";
        
    }
    // End profiling
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "x [m]\t y [mm]\n";
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << std::fixed << std::setprecision(3)
                  << x[i] << "\t" << y[i]*1000 << "\n";
    }
    std::cout << "\nN = " << N << "\n";
    double yAnal = solver.analyticMaxDeflection();
    std::cout << "\nAnalytic max deflection = " << yAnal*1000 << " mm\n";
    std::cout << "Numerical max deflection = " << (*std::max_element(y.begin(), y.end()))*1000 << " mm\n";

    // Výstup do CSV
    std::ofstream outFile(fileName);
    outFile << "x[m],y[m]\n";
    for (size_t i = 0; i < y.size(); ++i) {
        outFile << std::fixed << std::setprecision(6)
                << x[i] << "," << y[i] << "\n";
    }
    outFile.close();
    std::cout << "\nResults written to " << std::filesystem::current_path() << "/" << fileName << " \n";
    
    std::cout << "Simulation time: " << elapsed.count() << " s\n";

    
}