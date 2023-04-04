#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* data_file)
: 
{
}

// Condition initiale
double Function::foncf(const double x, const double y) const
{
  return 2*(y − y*y + x − x*x) ;
}

// Solution exacte
double Function::foncg(const double x, const double y) const
{
  return 0;
}

// Terme source
double Function::fonch(const double x, const double y) const
{
  return 0;
}


#define _FUNCTION_CPP
#endif
