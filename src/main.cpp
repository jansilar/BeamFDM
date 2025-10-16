#include "BeamSolver.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>

int main() {
    int N = 200;          // počet uzlů
    std::string fileName = "beam_results.csv";
    // Parametry nosníku
    double L = 2.0;           // délka [m]
    double E = 1.0e10;        // Young [Pa]
    double d = 0.10;          // průměr [m]
    double I = M_PI * std::pow(d, 4) / 64.0; // moment setrvačnosti [m^4]
    double q = 300.0;         // zatížení [N/m]

    BeamSolver solver(N, L, E, I, q);
    // Start profiling
    auto start = std::chrono::high_resolution_clock::now();
    // Run solver
    solver.solveStatic();
    // End profiling
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;


    const auto& y = solver.getY();
    const auto& x = solver.getX();

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