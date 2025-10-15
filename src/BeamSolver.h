#pragma once
#include <vector>
#include <cmath>
#include <iostream>

class BeamSolver {
public:
    BeamSolver(int N, double L, double E, double I, double q);

    // provede výpočet průhybu
    void solve();

    // vrátí výsledky
    const std::vector<double>& getY() const { return y; }
    const std::vector<double>& getX() const { return x; }

    // analytické srovnání pro konzolu
    double analyticMaxDeflection() const;

private:
    int N;              // počet uzlů
    double L;           // délka nosníku [m]
    double E;           // Youngův modul [Pa]
    double I;           // moment setrvačnosti [m^4]
    double q;           // rozložené zatížení [N/m]
    double h;           // krok sítě [m]
    std::vector<double> x; // uzly
    std::vector<double> y; // výsledek průhybu
};