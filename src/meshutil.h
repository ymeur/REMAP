#include <vector>
#include "elt.h"

void cptEltGeom(Elt& elt, const Coord &pole);
void cptAllEltsGeom(Elt *elt, int N, const Coord &pole);
void update_baryc(Elt *elt, int N);
void computeGradients(Elt **elt, int N);

