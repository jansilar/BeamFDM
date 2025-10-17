#include "DataWriter.h"

DataWriter::DataWriter(const std::string& filename) {
    outFile.open(filename);
    if (!outFile.is_open()) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
}

DataWriter::~DataWriter() {
    if (outFile.is_open()) {
        outFile.close();
    }
}

void DataWriter::writeHeader(const std::vector<double>& x) {
    outFile << "Time[s], ";
    for (size_t i = 0; i < x.size(); ++i) {
        outFile << "x="<<x[i];
        if (i < x.size() - 1)
            outFile << ", ";
    }
    outFile << "\n";
}

void DataWriter::writeStep(double t, const std::vector<double>& y) {
    outFile << t << ", ";
    for (size_t i = 0; i < y.size(); ++i) {
        outFile << y[i];
        if (i < y.size() - 1)
            outFile << ", ";
    }
    outFile << "\n";
}
