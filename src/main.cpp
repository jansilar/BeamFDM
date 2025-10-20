#include "BeamSolver.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>
#include <string>
#include <filesystem>
#include <chrono>
#include "DataWriter.h"


int main() {

    int N = 10;          // number of nodes
    std::string fileNameStatic = "beam_static_results.csv";
    std::string fileNameDynamic = "beam_dynamic_results.csv";
    // Beam parameters
    double L = 2.0;           // length [m]
    double E = 1.0e10;        // Young [Pa]
    double d = 0.10;          // diameter [m]
    double I = M_PI * std::pow(d, 4) / 64.0; // moment of inertia [m^4]
    double q = 300.0;         // external load [N/m]
    // Parameters for dynamic simulation
    double rho = 700.0;      // density [kg/m^3]
    double A = M_PI * std::pow(d, 2) / 4.0; // cross-section area [m^2]
    double omega1 = std::pow(1.875,2) * std::sqrt(E*I/(rho*A*std::pow(L,4))); // natural frequency of a clamped beam
    double ksi = 0.05;        // critical damping ratio (5%)
    double c = ksi*(2*omega1*rho*A);        // damping coefficient
    double dx = L / (N - 1);          // space step [m]
    double dt = 0.1 * std::pow(dx,2) / 2 * std::sqrt(rho*A/(E*I));        // time step [s] (stability condition)
    if (!(dt > 0.0)) {
        throw std::runtime_error("Computed dt is not positive.");
    }
    BeamSolver solver(N, L, E, I, rho, A, c);
    const auto& x = solver.getX();
    const auto& y = solver.getY();

    // Start profiling
    auto start = std::chrono::high_resolution_clock::now();

    //===================
    //static calculation:
    //===================
    std::cout << "======= Running static (steady-state) simulation... =======\n";
    // Solve static model
    solver.solveStatic(q);
    // Write results to console
    std::cout << "x [m]\t y [mm]\n";
    for (size_t i = 0; i < y.size(); ++i) {
        std::cout << std::fixed << std::setprecision(3)
                  << x[i] << "\t" << y[i]*1000 << "\n";
    }
    std::cout << "\nN = " << N << "\n";

    //Compare to the analytical solution
    double yAnal = solver.analyticMaxDeflection(q);
    std::cout << "\nAnalytic max deflection = " << yAnal*1000 << " mm\n";
    std::cout << "Numerical max deflection = " << (*std::max_element(y.begin(), y.end()))*1000 << " mm\n";

    // Write results to CSV
    std::ofstream outFile(fileNameStatic);
    outFile << "x[m],y[m]\n";
    for (size_t i = 0; i < y.size(); ++i) {
        outFile << std::fixed << std::setprecision(6)
                << x[i] << "," << y[i] << "\n";
    }
    outFile.close();
    std::string s = std::filesystem::current_path().string();
    std::cout << "\nStatic results written to " << std::filesystem::current_path().string() << "/" << fileNameStatic << " \n";
    std::cout << "======= Static simulation finished. =======\n\n\n";


    //===================
    //dynamic simulation:
    //===================
    std::cout << "======= Running dynamic simulation... =======\n";
    
    // Output to CSV - inicialization
    DataWriter dataWriter(fileNameDynamic);
    dataWriter.writeHeader(x);
    double finalTime = 1; // Total simulation time [s]
    double t = 0;
    solver.resetState(); // reset dynamic state to zeros
    dataWriter.writeStep(t, q, y); // write initial state
    double outputStep = 0.001; // data output time step  [s]
    int nInternalSteps = static_cast<int>(outputStep / dt);
    double qPulse = 0;
    // simulation loop
    while (t < finalTime){
        // external force actes between 0.03s and 0.06s
        qPulse = (0.03 < t && t < 0.06) ? q : 0;
        // perform multiple internal time steps
        solver.stepDynamic(qPulse, dt, nInternalSteps);
        t += dt*nInternalSteps;
        // write current state to CSV
        dataWriter.writeStep(t, qPulse, y);
    }
    std::cout << "Dynamic results written to " << std::filesystem::current_path() << "/" << fileNameDynamic << " \n";  
    std::cout << "======= Dynamic simulation finished. =======\n";
    // End profiling
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation took: " << elapsed.count() << " s\n";
    
}