#ifndef  __MAPPER_HPP__
#define __MAPPER_HPP__
#include "parallel_tree.hpp"

namespace sphereRemap {

enum verbosity
{
       SILENT = 0,
       PROGRESS = 1,
       DETAILS = 2
};

void cptOffsetsFromLengths(const int *lengths, int *offsets, int sz);

class Mapper
{
public:
       Mapper(int verbose = 2) : verbose(verbose), neighbourElements(NULL) {}
       ~Mapper();
       
       double buildSSTree(vector<Node>& srcMsh, vector<Node>& trgMsh)
       {
               sstree.build(srcMsh, trgMsh);
       }

       /** @param trgElts are the elements of the unstructured target grid
           Returns the timings for substeps: */
       vector<double> computeWeights(vector<Elt>& trgElts, vector<Elt>& srcElts, int interpOrder);

       /* where weights are returned after call to `computeWeights` */
       double *remapMatrix;
       int *srcAddress;
       int *srcRank;
       int *dstAddress;
       int nWeights;

private:
       /** @return number of weights (local to cpu) */
       int remap(Elt* elements, int nbElements, int order);

       void buildMeshTopology();
       void computeGrads();
       void computeIntersection(Elt* elements, int nbElements);

       int verbose;

       /** Holds adaptional leaf nodes as ghost cells for gradient computations (hold neighbour elements from other ranks).
           They will be inserted to the local trees but not saved in its leaf array */
       vector<Node> neighbourNodes; 

       int nbNeighbourElements;
       Elt* neighbourElements;

       CParallelTree sstree;
};

}
#endif
