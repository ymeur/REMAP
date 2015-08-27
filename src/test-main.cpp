#include <mpi.h>
#include <netcdf.h>
#include "errhandle.hpp"
#include "libmapper.hpp"

#include "node.hpp"

using namespace sphereRemap ;

int ncread(double *lon, double *lat, int nElt, const char* filename)
{
	int ncid, blonid, blatid, valid, neltid;
	size_t nEltCheck;
	exit_on_failure(nc_open(filename, NC_NOWRITE, &ncid), std::string("Cannot open netCDF file ") + filename);
	exit_on_failure(nc_inq_dimid(ncid, "elt", &neltid), std::string("No dimension elt in file ") + filename);
	nc_inq_dimlen(ncid, neltid, &nEltCheck);
	exit_on_failure(nElt != nEltCheck, std::string("Array sizes do not match!"));

	nc_inq_varid(ncid, "bounds_lon", &blonid);
	nc_inq_varid(ncid, "bounds_lat", &blatid);
	nc_get_var_double(ncid, blonid, lon);
	nc_get_var_double(ncid, blatid, lat);
	nc_close(ncid);
}

void compute_distribution(int nGlobalElts, int &start, int &nLocalElts)
{
	int mpiSize, mpiRank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
	start = 0;
	nLocalElts = 0;
	for (int i = 0; i <= mpiRank; i++)
	{
		start += nLocalElts;
		nLocalElts = nGlobalElts/mpiSize;
		if (i < nGlobalElts % mpiSize) nLocalElts++;
	}
}

int main()
{
	int interpOrder = 2;
/* low resolution
	char srcFile[] = "h14.nc";
	char dstFile[] = "r180x90.nc";
	double srcPole[] = {0, 0, 0};
	double dstPole[] = {0, 0, 1};
	int nSrcElt = 13661;
	int nDstElt = 16200;
*/
/* high resolution */
	char srcFile[] = "t740.nc";
	char dstFile[] = "r1440x720.nc";
	double srcPole[] = {0, 0, 0};
	double dstPole[] = {0, 0, 1};
	int nSrcElt = 741034;
	int nDstElt = 1036800;
         
	int nVert = 10;
	int nSrcLocal, nDstLocal, startSrc, startDst;
	int nWeight;

	MPI_Init(NULL, NULL);

	double *srcLon = (double *) malloc(nSrcElt*nVert*sizeof(double));	
	double *srcLat = (double *) malloc(nSrcElt*nVert*sizeof(double));	
	double *dstLon = (double *) malloc(nDstElt*nVert*sizeof(double));	
	double *dstLat = (double *) malloc(nDstElt*nVert*sizeof(double));	
	ncread(srcLon, srcLat, nSrcElt, srcFile);
	ncread(dstLon, dstLat, nDstElt, dstFile);

	compute_distribution(nSrcElt, startSrc, nSrcLocal);
	compute_distribution(nDstElt, startDst, nDstLocal);

	remap_get_num_weights(srcLon + startSrc*nVert, srcLat + startSrc*nVert, nVert, nSrcLocal, srcPole,
	                      dstLon + startDst*nVert, dstLat + startDst*nVert, nVert, nDstLocal, dstPole,
	                      interpOrder, &nWeight);
/*
	double *weights = (double *) malloc(nWeight*sizeof(double));	
	int    *srcIdx  = (int *)    malloc(nWeight*sizeof(int));	
	int    *srcRank = (int *)    malloc(nWeight*sizeof(int));	
	int    *dstIdx  = (int *)    malloc(nWeight*sizeof(int));	

	remap_get_weights(weights, srcIdx, dstIdx);

	free(srcLon); free(srcLat);
	free(dstLon); free(dstLat);
	free(srcIdx); free(dstIdx);
	free(srcRank);
	free(weights);

#ifdef DEBUG
	memory_report();
#endif
*/

	MPI_Finalize();
	return 0;
}
