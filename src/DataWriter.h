#pragma once
#include <fstream>
#include <vector>
#include <string>


/// @brief Object to write simulation data to a CSV file
class DataWriter {
public:
    /// @brief Constructor opens the file for writing
    /// @param filename Name of the CSV file
    DataWriter(const std::string& filename);
    ~DataWriter();

    /// @brief Writes the header row to the CSV file
    /// @param x Vector of x positions [m]
    void writeHeader(const std::vector<double>& x);

    /// @brief Writes a simulation step to the CSV file
    /// @param t Current time [s]
    /// @param q Current distributed load [N/m]
    /// @param y Current deflection values [m]
    void writeStep(double t, double q, const std::vector<double>& y);
private:
    std::ofstream outFile;
};