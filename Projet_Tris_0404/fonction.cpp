#include "fonction.h"
#include <cmath>

double f_cas1(double x, double y, double t, double Lx, double Ly)
{
    double res;
    
    res=2*(y-y*y+x-x*x);
    
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