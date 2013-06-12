/***************************************************************************
  ComBasic.h  -  Library to handle Basic Communication
  ------------------
 ***************************************************************************/

#ifndef ComBasic_H
#define ComBasic_H

#include <mpi.h>
#include "../utility/errors.h"
#include "../utility/debug.h"
#include "ComParser3D.h"


/** communicate particles along a direction **/
inline void communicateParticlesDIR(int buffer_size, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *b_right, double *b_left) {
  MPI_Status status;
  double *LEN = new double[3];
  LEN[0] = XLEN, LEN[1] = YLEN, LEN[2] = ZLEN;
  switch (DIR) {
    case 0:
      myrank = (int) floor((double) myrank / (ZLEN * YLEN));
      break;
    case 1:
      myrank = (int) floor((double) myrank / ZLEN);
      break;
  }
  if (myrank % 2 == 0) {
    // On the boundaries send e receive only if you have periodic condition: send to X-RIGHT
    if (right_neighbor != MPI_PROC_NULL) {
      if (LEN[DIR] > 1)
        MPI_Sendrecv_replace(&b_right[0], buffer_size, MPI_DOUBLE, right_neighbor, 1, right_neighbor, 1, MPI_COMM_WORLD, &status);
      else
        swapBuffer(buffer_size, b_left, b_right);
    }
  }
  else {
    // On the boundaries send e receive only if you have periodic condition: send to X-LEFT
    if (left_neighbor != MPI_PROC_NULL) {
      if (LEN[DIR] > 1)
        MPI_Sendrecv_replace(&b_left[0], buffer_size, MPI_DOUBLE, left_neighbor, 1, left_neighbor, 1, MPI_COMM_WORLD, &status);
      else
        swapBuffer(buffer_size, b_left, b_right);
    }
  }
  if (myrank % 2 == 1) {
    // On the boundaries send e receive only if you have periodic condition: send to X-RIGHT
    if (right_neighbor != MPI_PROC_NULL) {
      if (LEN[DIR] > 1)
        MPI_Sendrecv_replace(&b_right[0], buffer_size, MPI_DOUBLE, right_neighbor, 1, right_neighbor, 1, MPI_COMM_WORLD, &status);

    }
  }
  else {
    // On the boundaries send e receive only if you have periodic condition: send to X-LEFT
    if (left_neighbor != MPI_PROC_NULL) {
      if (LEN[DIR] > 1)
        MPI_Sendrecv_replace(&b_left[0], buffer_size, MPI_DOUBLE, left_neighbor, 1, left_neighbor, 1, MPI_COMM_WORLD, &status);
    }
  }
  delete[]LEN;
}
/** communicate ghost along a direction **/
inline void communicateGhostFace(int b_len, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightFace, double *ghostLeftFace) {

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
  if (rankF % 2 == 0 && right_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)  // 
    MPI_Sendrecv_replace(&ghostRightFace[0], b_len, MPI_DOUBLE, right_neighbor, 1, right_neighbor, 1, MPI_COMM_WORLD, &status);
  else if (rankF % 2 == 1 && left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)  // 
    MPI_Sendrecv_replace(&ghostLeftFace[0], b_len, MPI_DOUBLE, left_neighbor, 1, left_neighbor, 1, MPI_COMM_WORLD, &status);

  if (rankF % 2 == 1 && right_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)  // 
    MPI_Sendrecv_replace(&ghostRightFace[0], b_len, MPI_DOUBLE, right_neighbor, 1, right_neighbor, 1, MPI_COMM_WORLD, &status);
  else if (rankF % 2 == 0 && left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)  // 
    MPI_Sendrecv_replace(&ghostLeftFace[0], b_len, MPI_DOUBLE, left_neighbor, 1, left_neighbor, 1, MPI_COMM_WORLD, &status);
  // just swap the buffer if you have just a1 processor in 1 direction
  if (LEN[DIR] == 1 && right_neighbor != MPI_PROC_NULL && left_neighbor != MPI_PROC_NULL)
    swapBuffer(b_len, ghostLeftFace, ghostRightFace);


  delete[]LEN;
}

inline int my_MPI_Sendrecv(double *buf, int count,
              int dest, int source,
              MPI_Request *requests)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Datatype datatype = MPI_DOUBLE;
   int sendtag = 1;
   int recvtag = 1;
   //size_t numbytes = MPI_Sizeof(datatype)*count;
   //size_t numbytes = sizeof(double)*count;
   //void* sendbuffer = malloc(numbytes);
   //memcpy(sendbuffer, buf, numbytes);
   double* sendbuffer = &buf[count];
   for(int i=0;i<count;i++) sendbuffer[i]=buf[i];
   MPI_Isend(&sendbuffer[0], count, datatype,
     dest, sendtag, comm, &requests[0]);
   MPI_Irecv(buf, count, datatype,
     source, recvtag, comm, &requests[1]);
   MPI_Waitall(2,requests,MPI_STATUS_IGNORE);
   //free(sendbuffer);
}

inline int my_MPI_SendrecvNowait(double *buf, int count,
              int dest, int source,
              MPI_Request *requests)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Datatype datatype = MPI_DOUBLE;
   int sendtag = 1;
   int recvtag = 1;
   //size_t numbytes = MPI_Sizeof(datatype)*count;
   //size_t numbytes = sizeof(double)*count;
   //void* sendbuffer = malloc(numbytes);
   //memcpy(sendbuffer, buf, numbytes);
   double* sendbuffer = &buf[count];
   for(int i=0;i<count;i++) sendbuffer[i]=buf[i];
   MPI_Isend(&sendbuffer[0], count, datatype,
     dest, sendtag, comm, &requests[0]);
   MPI_Irecv(buf, count, datatype,
     source, recvtag, comm, &requests[1]);
   //MPI_Waitall(2,requests,MPI_STATUS_IGNORE);
   //free(sendbuffer);
}

/** communicate ghost along a direction **/
inline void communicateGhostNowait(int b_len,
  int DIR, double *ghostRightFace, double *ghostLeftFace,
  MPI_Request* requests,
  VirtualTopology3D* vct)
{
  MPI_Status status;
  const int myrank = vct->getCartesian_rank();
  int rankF;
  int const * const LEN = vct->getLEN();
  switch (DIR) {
    case 0:
      rankF = (int) floor((double) myrank / (LEN[2] * LEN[1]));
      break;
    case 1:
      rankF = (int) floor((double) myrank / LEN[2]);
      break;
    case 2:
      rankF = myrank;
      break;
  }
  int rght_neighbor=-1;
  int left_neighbor=-1;
  switch(DIR)
  {
    case 0:
      rght_neighbor = vct->getXright_neighbor_P();
      left_neighbor = vct->getXleft_neighbor_P();
      break;
    case 1:
      rght_neighbor = vct->getYright_neighbor_P();
      left_neighbor = vct->getYleft_neighbor_P();
      break;
    case 2:
      rght_neighbor = vct->getZright_neighbor_P();
      left_neighbor = vct->getZleft_neighbor_P();
      break;
    default:
      invalid_value_error(DIR);
  }
  if (rankF % 2 == 0) // order in the even-rank case:
  {
    if (rght_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)  // 
      my_MPI_SendrecvNowait(&ghostRightFace[0], b_len, rght_neighbor, rght_neighbor, &requests[0]);
    if (left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)  // 
      my_MPI_SendrecvNowait(&ghostLeftFace[0], b_len, left_neighbor, left_neighbor, &requests[2]);
  }
  else // order in the odd-rank case:
  {
    if (left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)
      my_MPI_SendrecvNowait(&ghostLeftFace[0], b_len, left_neighbor, left_neighbor, &requests[2]);
    if (rght_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)
      my_MPI_SendrecvNowait(&ghostRightFace[0], b_len, rght_neighbor, rght_neighbor, &requests[0]);
  }

  // just swap the buffer if you have just a1 processor in 1 direction
  if (LEN[DIR] == 1 && rght_neighbor != MPI_PROC_NULL && left_neighbor != MPI_PROC_NULL)
    swapBuffer(b_len, ghostLeftFace, ghostRightFace);
}

/** non-waiting communicate ghost along a direction **/
inline void communicateGhostNowait_broken(int b_len, int DIR,
  double *rghtBuffer, double *leftBuffer,
  MPI_Request* requests,
  VirtualTopology3D* vct)
{
  //LEN[0] = XLEN, LEN[1] = YLEN, LEN[2] = ZLEN;
  int const * const LEN = vct->getLEN();

  // get the neighbors in the relevant direction.
  // should implement
  // VirtualTopology3D::getLeft_neighbor_P(int DIR) and
  // VirtualTopology3D::getRight_neighbor_P(int DIR).
  //
  int rght_neighbor=-1;
  int left_neighbor=-1;
  switch(DIR)
  {
    case 0:
      rght_neighbor = vct->getXright_neighbor_P();
      left_neighbor = vct->getXleft_neighbor_P();
      break;
    case 1:
      rght_neighbor = vct->getYright_neighbor_P();
      left_neighbor = vct->getYleft_neighbor_P();
      break;
    case 2:
      rght_neighbor = vct->getZright_neighbor_P();
      left_neighbor = vct->getZleft_neighbor_P();
      break;
    default:
      invalid_value_error(DIR);
  }
  if (rght_neighbor != MPI_PROC_NULL && LEN[DIR] > 1) {
    for(int i=0;i<b_len;i++)
    {
      rghtBuffer[b_len+i]=rghtBuffer[i];
    }
    MPI_Isend(&rghtBuffer[b_len], b_len, MPI_DOUBLE, rght_neighbor, 1, MPI_COMM_WORLD, &requests[0]);
    MPI_Irecv(&rghtBuffer[0], b_len, MPI_DOUBLE, rght_neighbor, 1, MPI_COMM_WORLD, &requests[1]);
  }

  if (left_neighbor != MPI_PROC_NULL && LEN[DIR] > 1)
  {
    for(int i=0;i<b_len;i++)
    {
      leftBuffer[b_len+i]=leftBuffer[i];
    }
    MPI_Isend(&leftBuffer[b_len], b_len, MPI_DOUBLE, left_neighbor, 1, MPI_COMM_WORLD, &requests[2]);
    MPI_Irecv(&leftBuffer[0], b_len, MPI_DOUBLE, left_neighbor, 1, MPI_COMM_WORLD, &requests[3]);
  }

  // just swap the buffer if you have just a single processor in the direction
  if (LEN[DIR] == 1 && rght_neighbor != MPI_PROC_NULL && left_neighbor != MPI_PROC_NULL)
    swapBuffer(b_len, leftBuffer, rghtBuffer);
}

inline void waitForRequests(MPI_Request* requests, int* _request_start, int request_idx)
{
  int& request_start = *_request_start;
  // make sure the most recent communication is completed
  //
  MPI_Waitall(4*(request_idx-request_start),
    &requests[4*request_start],
    MPI_STATUS_IGNORE);
  request_start = request_idx;
}

#endif
