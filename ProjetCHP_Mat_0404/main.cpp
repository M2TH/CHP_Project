// Code principale du projet de C++

//compilation : 



// #include <mpi.h>
#include <stdio.h>  //pour le printf
#include <iostream>


int main(int argc, char ** argv){

    //using namespace std;

     // Ouverture du fichier de paramètre
    std::ifstream param_file("parameters.txt");
    if (!param_file.is_open())
    {
        std::cerr << "Error: Could not open parameter file.\n";
        return 1;
    }

    // Read the parameters from the file
    int Nx, Ny;
    double Lx, Ly, D, dt;
    std::string line;

    while (std::getline(param_file, line))
    {
        if (line[0] == '#')
        {
            // Skip comment lines starting with #
            continue;
        }
        else if (line.find("Nx") != std::string::npos)
        {
            // Read Nx parameter
            std::istringstream iss(line.substr(line.find("=") + 1));
            iss >> Nx;
        }
        else if (line.find("Ny") != std::string::npos)
        {
            // Read Ny parameter
            std::istringstream iss(line.substr(line.find("=") + 1));
            iss >> Ny;
        }
        else if (line.find("Lx") != std::string::npos)
        {
            // Read Lx parameter
            std::istringstream iss(line.substr(line.find("=") + 1));
            iss >> Lx;
        }
        else if (line.find("Ly") != std::string::npos)
        {
            // Read Ly parameter
            std::istringstream iss(line.substr(line.find("=") + 1));
            iss >> Ly;
        }
        else if (line.find("D") != std::string::npos)
        {
            // Read D parameter
            std::istringstream iss(line.substr(line.find("=") + 1));
            iss >> D;
        }
        else if (line.find("dt") != std::string::npos)
        {
            // Read dt parameter
            std::istringstream iss(line.substr(line.find("=") + 1));
            iss >> dt;
        }
    }

    // Close the parameter file
    param_file.close();

    // Print the read parameters to console (optional)
    std::cout << "Nx = " << Nx << std::endl;
    std::cout << "Ny = " << Ny << std::endl;
    std::cout << "Lx = " << Lx << std::endl;
    std::cout << "Ly = " << Ly << std::endl;
    std::cout << "D = " << D << std::endl;
    std::cout << "dt = " << dt << std::endl;


    return 0;
}