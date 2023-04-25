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
#include "gradConj.h"
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
    



    int cas, maxiter, n;
    double alpha,beta,gamma,dx,dy,tol;
    double x,y,t;
 

    cas=2; //cas d'etude
    tol=0.01;
    maxiter=1000;

    
    
    n=Nx*Ny; // dim du maillage
    dx=Lx/(Nx+1); // pas selon x
    dy=Ly/(Ny+1); // pas selon y
    
    alpha=1+(2/(dx*dx))+(2/(dy*dy));
    beta=-1/(dx*dx);
    gamma=-1/(dy*dy);

    vector<double> X(n),res_vec(n),condInit(n),b(n);
   
    
    for(int i=0;i<n;i++)
    {
        X[i]=0;
    }

    X[4]=1;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // Test du produit matrice vecteur
    /////////////////////////////////////////////////////////////////////////////////////////////

    // res_vec=matvec(alpha,beta,gamma,Nx,Ny,X);
    // cout<<"alpha="<<alpha<<" beta="<<beta<<" gamma="<<gamma<<endl;
    // cout<<"Voici le résultat :                      avec X:"<<endl;
    // for(int i=0;i<n;i++)
    // {
    //     cout<<res_vec[i]<<"                         "<<X[i]<<endl;
    // }



    // for(int i=0;i<n;i++)
    // {
    //     X[i]=0;
    // }

    //X[0]=1;

    t=0;
    x=dx;
    y=dy;
    cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;
  
    res_vec=F_b(Lx,Ly,beta,gamma,dt,t,D,cas,Nx,Ny,X);
    // cout<<"dx="<<dx<<" , dy="<<dy<<endl;
    // cout<<dt*f(x,y,t,Lx,Ly,cas)<<endl;
    // cout<<"Voici le second membre :                      avec X:"<<endl;
    // for(int i=0;i<n;i++)
    // {
    //    cout<<res_vec[i]<<"                         "<<X[i]<<endl;
    // }

    cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;

    /////////////////////////////////////////////////////////////////////////////////////////////
    // Coeur du programme : appel du gradient conjugué
    /////////////////////////////////////////////////////////////////////////////////////////////



    // Construction d'un vecteur condition initial
    for(int i=0;i<n;i++)
    {
        condInit[i]=1;
    }

    cout<<"//////////////////////////////////////////////////////////////////////////"<<endl;

    X = gradientConj( condInit , tol, maxiter, Nx, Ny, alpha, beta, gamma,res_vec);

    

    // cout<<"Voici la condition initial :                      avec X:"<<endl;
    // for(int i=0;i<n;i++)
    // {
    //     cout<<condInit[i]<<"                         "<<X[i]<<endl;
    // }

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    /// mettre dans un fichier
    ////////////////////////////////////////////////////////////////////////////////////////////////////

    
    dt=1;
    cas=2;
    
    for(int k=0;k<10;k++)
    {
        t=k*dt;
        b=F_b(Lx,Ly,beta,gamma,dt,t,D,cas,Nx,Ny,X);
        X = gradientConj( X , tol, maxiter, Nx, Ny, alpha, beta, gamma,res_vec);

        ofstream mon_flux; // Contruit un objet "ofstream"
        string name_file("sol.10.dat"); // Le nom de mon fichier
        mon_flux.open(name_file, ios::out); // Ouvre un fichier appelé name_file


        /////////////////////////////////////////////////////////////////
        /// cas pour i différent de 0 et Ny-1
        /////////////////////////////////////////////////////////////////
        for (int i=0;i<Ny;i++) 
        {
            y=(i+1)*dy;
            for(int j=0;j<Nx;j++) // cas pour j différent de 0 et Nx-1
            {
                mon_flux << (1+j)*dx << " " << y << " " << X[j+i*Nx] << " " << endl;
            }
        }
        
        mon_flux.close(); // Ferme le fichier
    }

    return 0;
}