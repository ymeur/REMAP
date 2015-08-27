#ifndef  __TREE_H__
#define __TREE_H__
#include <iostream>
#include "triple.hpp"

namespace sphereRemap {

struct CGrid
{
	Coord pole;
	int numElts;
};

Coord readPole(std::istream&);

extern CGrid srcGrid;
extern CGrid tgtGrid;

}

#endif
