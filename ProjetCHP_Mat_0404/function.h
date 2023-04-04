#ifndef _FUNCTION_H

// #include "DataFile.h"

class Function {
private:
  // Quelques variables privées utiles pour
  // construire la condition initiale et la solution exacte (sigma)
  
  //const double _xmin, _ymin;
  // const double _xmax, _ymax;

public:
  // Constructeur
  Function(DataFile* data_file);

  // fonction f
  double foncf(const double x, const double y) const;

  // fonction g
  double foncg(const double x, const double y) const;

  //fonction h
  double fonch(const double x, const double y) const;

};

#define _FUNCTION_H
#endif