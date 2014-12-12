#include <mpi.h>
#include "ComBasic3D.h"
#include "ComParser3D.h"
#include "ipic_defs.h"
#include "TimeTasks.h"
#include "debug.h"
#include <math.h>

/** communicate ghost along a direction **/
void communicateGhostFace(int b_len, int myrank,
  int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN,
  double *ghostRightFace, double *ghostLeftFace)
{
  timeTasks_set_communicating();
  MPI_Status status;
  double *LEN = new double[3];
  int rankF;
  LEN[0] = XLEN, LEN[1] = YLEN, LEN[2] = ZLEN;
  switch (DIR) {
    case 0:
      rankF = (int) floor((double) myrank / (ZLEN * YLEN));
      break;
    case 1:
      rankF = (int) floor((double) myrank / ZLEN);
      break;
    case 2:
      rankF = myrank;
      break;
  }
  if (rankF % 2 == 0 && right_neighbor != MPI_PROC_NULL && LEN[DIR] > 1){
    ipic_MPI_Sendrecv_replace(&ghostRightFace[0], b_len, MPI_DOUBLE,
      right_neighbor, 1, right_neighbor, 1, MPI_COMM_WORLD, &status);
  }else if (rankF % 2 == 1 && left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1){
    ipic_MPI_Sendrecv_replace(&ghostLeftFace[0], b_len, MPI_DOUBLE,
      left_neighbor, 1, left_neighbor, 1, MPI_COMM_WORLD, &status);
  }

  if (rankF % 2 == 1 && right_neighbor != MPI_PROC_NULL && LEN[DIR] > 1){
    ipic_MPI_Sendrecv_replace(&ghostRightFace[0], b_len, MPI_DOUBLE,
      right_neighbor, 1, right_neighbor, 1, MPI_COMM_WORLD, &status);
  }else if (rankF % 2 == 0 && left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1){
    ipic_MPI_Sendrecv_replace(&ghostLeftFace[0], b_len, MPI_DOUBLE,
      left_neighbor, 1, left_neighbor, 1, MPI_COMM_WORLD, &status);
  }
  // just swap the buffer if you have just a1 processor in 1 direction
  if (LEN[DIR] == 1 && right_neighbor != MPI_PROC_NULL && left_neighbor != MPI_PROC_NULL)
    swapBuffer(b_len, ghostLeftFace, ghostRightFace);


  delete[]LEN;
}
