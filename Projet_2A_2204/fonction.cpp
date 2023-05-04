#include "fonction.h"
#include <cmath>
#include <vector>
#include <iostream>

//////////////////////////////////////////////////////////////////////////
//// fonction f,g et h en fonction des différents cas
//////////////////////////////////////////////////////////////////////////

double f(double x, double y, double t, double Lx, double Ly, int cas)
{
    // fonction f pour chacun des cas
    double res;
    if(cas==1)
    {
        res=2*(y-y*y+x-x*x);
    }
    else if(cas==2)
    {
        res=sin(x)+cos(y);
    }
    else if(cas==3)
    {
        res=exp(-pow(x-(Lx/2),2))*exp(-pow(y-(Ly/2),2))*cos(M_PI*t/2);
    }
    return res;
}

double g(double x, double y, double Lx, double Ly, int cas)
{
    // fonction g pour chacun des cas
    double res;
    if(cas==1)
    {
        res=0.0 ;
    }
    else if(cas==2)
    {
        res=sin(x)+cos(y);
    }
    else if(cas==3)
    {
        res=0.0 ;
    }
    return res;
}

double h(double x, double y, double Lx, double Ly, int cas)
{
    // fonction h pour chacun des cas
    double res;
    if(cas==1)
    {
        res=0.0 ;
    }
    else if(cas==2)
    {
        res=sin(x)+cos(y);
    }
    else if(cas==3)
    {
        res=1.0 ;
    }
    return res;
}

vector<double> F_b(double Lx,double Ly, double beta, double gamma, double dt, double t, double D, int cas, int Nx, int Ny, vector<double> U)
{
    // fonction pour avoir le vecteur b du second membre 
    vector <double> F(Nx*Ny);
    double x,y,dx,dy;

    dx=Lx/(Nx+1);
    dy=Ly/(Ny+1);

    x=dx;
    y=dy;
    
    ////////////////////////////////////////////////////////////////
    // boucle de définition du bord du bas (y=0 et x entre 1 et Nx) ------------------------------------------------------------------------
    ////////////////////////////////////////////////////////////////

    F[0]=U[0]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(0,y,Lx,Ly,cas)-dt*D*gamma*g(x,0,Lx,Ly,cas); // cas j=0 (j boucle sur x) = coin en bas à gauche
    
    for(int j=1;j<Nx-1;j++) // segment du bord du bas (sauf les deux coins) cas pour j différent de 0 et Nx-1
    {
        x=x+dx;
        F[j]=U[j]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*gamma*g(x,0,Lx,Ly,cas);
        
    }
    
    //cas j=Nx-1 : bord en bas à droite
    x=x+dx;
    F[Nx-1]=U[Nx-1]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(x+dx,y,Lx,Ly,cas)-dt*D*gamma*g(x,0,Lx,Ly,cas);
    


    /////////////////////////////////////////////////////////////////
    /// Serpent pour le bloc centrale et les conditions aux limites sur les bords gauches et droits: cas pour i différent de 0 et Ny-1
    /////////////////////////////////////////////////////////////////
    for (int i=1;i<Ny-1;i++) 
    {
        // cas j=0 : bord gauche (x=0)(sauf coin en bas et en haut à gauche)
        int df; // indice du premier F du bloc
        df=i*Nx;
        x=dx;
        y=y+dy;
        F[df]=U[df]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(0,y,Lx,Ly,cas);
        //cout<<"h="<<h(0,y,Lx,Ly,cas)<<"   F="<<F[df]<<"   beta="<<beta<<"  f="<<f(x,y,t,Lx,Ly,cas)<<endl;

        // centre du domaine (qd il n'y a pas d'influence des conditions de bords)
        for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
        {
            x=x+dx;
            F[df+j]=U[df+j]+dt*f(x,y,t,Lx,Ly,cas);
            //cout<<"df+j="<< df+j <<"   F="<<F[df+j] << endl;
        }

        // cas j=Nx-1 : bord droit (x=Lx)(sauf coin en bas et en haut à droite)
        x=x+dx;
        F[df+Nx-1]=U[df+Nx-1]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(Lx,y,Lx,Ly,cas); //F[df+Nx-1]=U[df+Nx-1]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(x+dx,y,Lx,Ly,cas);
        // cout<<"y="<< y <<"  df+Nx-1="<< df+Nx-1 <<"   F="<<F[df+Nx-1] << endl;
        // cout<<"h="<<h(x+dx,y,Lx,Ly,cas)<<"   F="<<F[df+Nx-1]<<"   beta="<<beta<<"  f="<<f(x,y,t,Lx,Ly,cas)<<endl;
    }

    /////////////////////////////////////////////////////////////////
    /// cas i=Ny-1 : bord du haut
    /////////////////////////////////////////////////////////////////
    
    // cas j=0
    int df;
    df=Nx*(Ny-1);
    x=dx;
    y=y+dy;
    
    F[df]=U[df]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(0,y,Lx,Ly,cas)-dt*D*gamma*g(x,y+dy,Lx,Ly,cas);

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
    {
        x=x+dx;
        F[df+j]=U[df+j]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*gamma*g(x,y+dy,Lx,Ly,cas);
    }
    
    // cas j=Nx-1
    x=x+dx;
    F[df+Nx-1]=U[df+Nx-1]+dt*f(x,y,t,Lx,Ly,cas)-dt*D*beta*h(Lx,y,Lx,Ly,cas)-dt*D*gamma*g(x,Ly,Lx,Ly,cas);
    
    return F;
}


//////////////////////////////////////////////////////////////////////////
//// fonction multiplication matrice A* vecteur U
//////////////////////////////////////////////////////////////////////////

vector<double> matvec(double alpha, double beta, double gamma, int Nx, int Ny, vector<double> U)
{
    vector<double> res(Nx*Ny);
    // alpha : diagonale central
    // beta : surdiagonale
    // gamma : diagonale externe
    /////////////////////////////////////////////////////////////////
    /// cas i=0 (i boucle sur les blocs)
    /////////////////////////////////////////////////////////////////

    res[0]=alpha*U[0]+beta*U[1]+gamma*U[Nx]; // cas j=0 (j boucle dans les blocs)
    
    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx-1
    {
        res[j]=beta*U[j-1]+alpha*U[j]+beta*U[j+1]+gamma*U[j+Nx];
    }
    
    res[Nx-1]=beta*U[Nx-2]+alpha*U[Nx-1]+gamma*U[Nx+Nx-1]; // cas j=Nx-1

    /////////////////////////////////////////////////////////////////
    /// cas pour i différent de 0 et Ny-1
    /////////////////////////////////////////////////////////////////
    for (int i=1;i<Ny-1;i++) 
    {
        int dy; // indice du alpha en haut à gauche du bloc
        dy=i*Nx;
        res[dy]=gamma*U[dy-Nx]+alpha*U[dy]+beta*U[dy+1]+gamma*U[dy+Nx]; // cas j=0

        for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
        {
            res[j+dy]=gamma*U[dy-Nx+j]+beta*U[dy+j-1]+alpha*U[dy+j]+beta*U[dy+j+1]+gamma*U[dy+j+Nx];
        }

        res[dy+Nx-1]=gamma*U[dy-1]+beta*U[dy+Nx-2]+alpha*U[dy+Nx-1]+gamma*U[dy+Nx]; // cas j=Nx-1
    }
    
    /////////////////////////////////////////////////////////////////
    /// cas i=Ny-1
    /////////////////////////////////////////////////////////////////
    int dy;
    dy=Nx*(Ny-1);

    res[dy]=gamma*U[dy-Nx]+alpha*U[dy]+beta*U[dy+1]; // cas j=0

    for(int j=1;j<Nx-1;j++) // cas pour j différent de 0 et Nx
    {
        res[dy+j]=gamma*U[dy+j-Nx]+beta*U[dy+j-1]+alpha*U[dy+j]+beta*U[dy+j+1];
    }

    res[dy+Nx-1]=gamma*U[dy-1]+beta*U[dy+Nx-2]+alpha*U[dy+Nx-1]; // cas j=Nx-1
    
    return res;
}

void charge(int me, int n, int np, int *iBeg, int *iEnd)
{
    int divi, rest;

    rest=n%np;
    divi=n/np;
    if(me<rest)
    {
        *iBeg=me*divi+me;
        *iEnd=(me+1)*divi+me;
    }
    else
    {
        *iBeg=me*divi+rest;
        *iEnd=(me+1)*divi-1+rest;
    }
}