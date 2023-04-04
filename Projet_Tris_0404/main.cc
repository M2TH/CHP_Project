#include <iostream>
#include <stdio.h>
#include <vector>
#include "fonction.cpp"

using namespace std;

int main()
{
    int nx,ny;
    double alpha,beta,gamma,D,Lx,Ly,deltax,deltay,deltat;
    double x,y,t;
    double res;
    
    x=2;
    y=1;
    t=1;
    Lx=1;
    Ly=1;
    nx=3;
    ny=4;
    deltax=Lx/nx;
    deltay=Ly/ny;
    deltat=1;
    alpha=(2/(deltax*deltax))*(2/(deltay*deltay));
    beta=-1/(deltax*deltax);
    gamma=-1/(deltay*deltay);

    vector<double> X;
    res=f_cas1(x,y,t,Lx,Ly);

    cout<<res<<endl;






    return 0;
}