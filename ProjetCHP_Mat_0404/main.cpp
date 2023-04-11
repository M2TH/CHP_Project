// Code principale du projet de C++

//compilation : 
// cf makefile


// #include <mpi.h>
#include <stdio.h>  //pour le printf
#include <iostream>
#include <fstream>


int main(int argc, char ** argv){

    using namespace std;

    // Read the parameters from the file
    int Nx, Ny;
    double Lx, Ly, D, dt;

    //Ouvre un fichier
    ifstream flux("parameter.txt");

    // Premiere facon de lire un fichier (ligne par ligne)
    string ligne;
    
    //Nx
    getline(flux, ligne);
    flux >> Nx ; 
    // getline(flux, value);
    // value = stoi(ligne);
    cout << ligne << endl << Nx << endl;

    





    // Close the parameter file
    flux.close();

    return 0;
}