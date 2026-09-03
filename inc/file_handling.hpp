// file_handling.hpp
#ifndef FILE_HANDLING_HPP
#define FILE_HANDLING_HPP

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

inline void write_csv(const std::string& filename,
                      const std::vector<double>& x,
                      const std::vector<double>& u,
                      const std::string& x_name,
                      const std::string& u_name)
{
    if (x.size() != u.size())
        throw std::runtime_error("x and u must have same size");

    std::ofstream file(filename);

    if (!file)
        throw std::runtime_error("couldn't open file: " + filename);

    file << std::setprecision(15) << std::scientific;
    file << x_name << "," << u_name << std::endl;

    for (size_t i = 0; i < x.size(); ++i)
        file << x[i] << "," << u[i] << "\n";
}

inline void read_csv(const std::string& filename,
                     std::vector<double>& x,
                     std::vector<double>& u)
{
    std::ifstream file(filename);

    if (!file)
        throw std::runtime_error("Couldn't open file: " + filename);

    x.clear();
    u.clear();

    std::string line;

    // skip header
    std::getline(file, line);

    double xi, ui;
    char sep;
    while (file >> xi >> sep >> ui)
    {
        if (sep != ',')
            throw std::runtime_error("invalid csv format");

        x.push_back(xi);
        u.push_back(ui);
    }
}

inline void write_line_to_file(const std::string& filename,
                              const std::string& line)
{
    std::ofstream file(filename, std::ios::app);
    if(!file.is_open())
        throw std::runtime_error("failed to open file: " + filename);

    file << line << "\n";
}

inline void print_csv(const std::string& filename)
{
    std::vector<double> x, u;
    read_csv(filename, x, u);

    std::cout << "Printing CSV file: " << filename << std::endl;

    std::cout << std::setprecision(15) << std::scientific;
    std::cout << "x,u\n";
    for (size_t i = 0; i < x.size(); ++i)
        std::cout << x[i] << "," << u[i] << "\n";
}
#endif // FILE_HANDLING_HPP
