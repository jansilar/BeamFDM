#pragma once
#include <vector>
#include <cmath>
#include <iostream>

class BeamSolver {
public:
    BeamSolver(int N, double L, double E, double I, double rho, double A, double c);

    // provede výpočet průhybu
    void solveStatic(double q);
    void stepDynamic(double q, double dt, int steps);

    // vrátí výsledky
    const std::vector<double>& getY() const { return y; }
    const std::vector<double>& getX() const { return x; }

    // analytické srovnání pro konzolu
    double analyticMaxDeflection(double q) const;

private:
    int N;              // počet uzlů
    double L;           // délka nosníku [m]
    double E;           // Youngův modul [Pa]
    double I;           // moment setrvačnosti [m^4]
    double dx;           // krok sítě [m]
    double rho;         // hustota materiálu [kg/m^3] - pro dynamiku
    double A;           // průřezová plocha [m^2] - pro dynamiku
    double c;           // tlumení - pro dynamiku
    std::vector<double> x; // uzly
    std::vector<double> y; // průhyb
    std::vector<double> yPrev; // průhyb v minulém kroce (pro dynamiku)
    std::vector<double> yNew; // průhyb v příštím kroce (pro dynamiku)
};