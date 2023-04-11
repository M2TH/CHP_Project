#include <iostream>
#include <stdio.h>
#include <vector>
#include <cmath>
#include "fonction.h"
 


int main()
{
    int Nx,Ny,n,cas;
    double alpha,beta,gamma,D,Lx,Ly,deltax,deltay,deltat;
    double x,y,t;
    double res;
    
    x=2; 
    y=1;
    t=1;
    Lx=1;
    Ly=1;
    Nx=4;
    Ny=3;
    n=Nx*Ny;
    cas=1;
    deltax=Lx/Nx;
    deltay=Ly/Ny;
    deltat=1;
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

