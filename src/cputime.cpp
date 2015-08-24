#include <mpi.h>

double cputime()
{
	return MPI_Wtime();
}
