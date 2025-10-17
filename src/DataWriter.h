#pragma once
#include <fstream>
#include <vector>
#include <string>

class DataWriter {
public:
    DataWriter(const std::string& filename);
    ~DataWriter();
    void writeHeader(const std::vector<double>& x);
    void writeStep(double t, double q, const std::vector<double>& y);
private:
    std::ofstream outFile;
};