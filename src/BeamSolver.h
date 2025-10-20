#pragma once
#include <vector>
#include <cmath>
#include <iostream>

/// @brief Class to solve static and dynamic beam deflection problems using FDM
class BeamSolver {
public:
    /// @brief Constructor
    /// @param N Number of nodes
    /// @param L Length of the beam [m] 
    /// @param E Young's modulus [Pa]
    /// @param I Moment of inertia [m^4]
    /// @param rho Density [kg/m^3]
    /// @param A Cross-sectional area [m^2]
    /// @param c Damping coefficient
    BeamSolver(int N, double L, double E, double I, double rho, double A, double c);

    /// @brief static (steady-state) deflection
    /// @param q Distributed load [N/m]
    void solveStatic(double q);

    /// @brief dynamic simulation steps
    /// @param q Distributed load [N/m]
    /// @param dt Time step [s]
    /// @param steps Number of time steps to perform
    void stepDynamic(double q, double dt, int steps);

    /// @brief Reset state to initial conditions (all zeros)
    void resetState();

    /// @brief Get results
    const std::vector<double>& getY() const { return y; }

    /// @brief  Get x coordinate values
    const std::vector<double>& getX() const { return x; }

    /// @brief solution of the static problem
    /// @param q Distributed load [N/m]
    double analyticMaxDeflection(double q) const;

private:
    int N;              // number of nodes
    double L;           // Beam length [m]
    double E;           // Young's modulus [Pa]
    double I;           // Moment of inertia [m^4]
    double dx;           // Space step [m]
    double rho;         // Density [kg/m^3] - (dynamics only)
    double A;           // Crosse [m^2] - pro dynamiku
    double c;           // tlumení - pro dynamiku
    std::vector<double> x; // x positions (nodes) [m]
    std::vector<double> y; // deflection values [m]
    std::vector<double> yPrev; // deflection in previou step (dynamics only)
    std::vector<double> yNew; // deflection in next step (dynamics only)
};