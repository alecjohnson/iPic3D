/***************************************************************************
  ComBasic.h  -  Library to handle Basic Communication
  ------------------
 ***************************************************************************/

#ifndef ComBasic_H
#define ComBasic_H

/** communicate ghost along a direction **/
void communicateGhostFace(int b_len, int myrank, int right_neighbor, int left_neighbor, int DIR, int XLEN, int YLEN, int ZLEN, double *ghostRightFace, double *ghostLeftFace);

#endif
