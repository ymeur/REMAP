#include <mpi.h>

namespace sphereRemap {
  
double cputime()
{
	return MPI_Wtime();
}

}
