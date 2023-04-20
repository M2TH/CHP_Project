// Code principale du projet de C++

//compilation : 
// cf makefile


// #include <mpi.h>
#include <stdio.h>  //pour le printf
#include <iostream>
#include <sstream> 
#include <fstream>
#include <vector>
#include <cmath>
#include "fonction.h"
#include <string>


int main(int argc, char ** argv){

    using namespace std;
    


    /////////////////////////////////////////
    // Lecture des parametres dans le fichier
    /////////////////////////////////////////

    int Nx, Ny;
    double Lx, Ly, D, dt;

    //Ouvre un fichier
    ifstream flux("parameter.txt");

    string ligne, valeur;
    
    //Nx
    getline(flux, ligne); //lecture de la chaine de caractères de la première ligne, stocké dans la chaine de caractère "ligne"
    getline(flux, valeur); //lecture de la chaine de caractères valeur
    istringstream iss(valeur); //copie de la chaine de caratères "valeur" dans l'objet iss de type std::istringstream
    iss >> Nx ;     // lecture de la valeur numérique à partir de iss
    cout << ligne << endl ; // affichage du paramètre et de sa valeur
    cout << Nx << endl;
    iss.clear(); // nettoyage de iss afin de pouvoir lui donner une autre valeur

    // Ny 
    getline(flux, ligne);
    getline(flux, valeur);
    iss.str(valeur); 
    iss >> Ny ;
    cout << ligne << endl ;
    cout << Ny << endl;
    iss.clear();

    // Lx 
    getline(flux, ligne);
    getline(flux, valeur);
    iss.str(valeur); 
    iss >> Lx ;
    cout << ligne << endl ;
    cout << Lx << endl;
    iss.clear();

    // Ly 
    getline(flux, ligne);
    getline(flux, valeur);
    iss.str(valeur); 
    iss >> Ly ;
    cout << ligne << endl ;
    cout << Ly << endl;
    iss.clear();

    // D 
    getline(flux, ligne);
    getline(flux, valeur);
    iss.str(valeur); 
    iss >> D ;
    cout << ligne << endl ;
    cout << D << endl;
    iss.clear();

    //dt 
    getline(flux, ligne);
    getline(flux, valeur);
    iss.str(valeur); 
    iss >> dt ;
    cout << ligne << endl ;
    cout << dt << endl;
    iss.clear();

    // Fermeture du fichier de parametre
    flux.close();
    



    int cas, n;
    double alpha,beta,gamma,deltax,deltay,deltat;
    double x,y,t;
    double res; // Qu'est ce que res ? (Inutile dans le main je crois)

    /*
    Lx=1;
    Ly=1;
    Nx=4;
    Ny=3;
    D = 1 ; */
    deltat=dt ; 

    cas=1; //cas d'etude

    x=2; // A quoi servent x, y et t  ? est ce la position en x, y et le temps courant respectivement? 
    y=1;
    t=1;
    
    n=Nx*Ny; // dim du maillage
    deltax=Lx/Nx; // pas selon x
    deltay=Ly/Ny; // pas selon y
    
    alpha=1+(2/(deltax*deltax))+(2/(deltay*deltay));
    beta=-1/(deltax*deltax);
    gamma=-1/(deltay*deltay);

    vector<double> X(n),res_vec(n);
   
    
    for(int i=0;i<n;i++)
    {
        X[i]=0;
    }

    X[4]=1;

    res_vec=matvec(alpha,beta,gamma,Nx,Ny,X);
    cout<<"alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<endl;
    cout<<"Voici le résultat :"<<endl;
    for(int i=0;i<n;i++)
    {
        cout<<res_vec[i]<<endl;
    }



    return 0;
}