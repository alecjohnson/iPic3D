
#include <mpi.h>
#include "ComParser3D.h"

/** swap the buffer */
void swapBuffer(int buffer_size, double *b_left, double *b_rght) {
  double *temp = new double[buffer_size];
  for (int i = 0; i < buffer_size; i++) {
    temp[i] = b_left[i];
    b_left[i] = b_rght[i];
    b_rght[i] = temp[i];
  }
  delete[]temp;
}
/** swap the buffer */
void swapBuffer(double *b_left, double *b_rght) {
  double temp = *b_left;
  *b_left = *b_rght;
  *b_rght = temp;

}
/** swap ghost cells */
void swapGhostFace(int n1, int n2, double **ghostFaceLeft, double **ghostFacerght) {
  double **temp = newArr2(double, n1, n2);
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      temp[i][j] = ghostFaceLeft[i][j];
      ghostFaceLeft[i][j] = ghostFacerght[i][j];
      ghostFacerght[i][j] = temp[i][j];
    }
  }
  delArr2(temp, n1);
}
/** prepare ghost cells on 6 faces for communication */
void makeNodeFace(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtFace, double *ghostXleftFace,
  double *ghostYrghtFace, double *ghostYleftFace,
  double *ghostZrghtFace, double *ghostZleftFace)
{
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // XFACES 
  int counter = 0;
  //#pragma omp for collapse(2)
  for (int j = 1; j <= nyr; j++)
  for (int k = 1; k <= nzr; k++)
  {
      const int count = (j-1)*nzr+k-1;
      assert_eq(count,counter);
      ghostXleftFace[count] = vector[2][j][k];
      ghostXrghtFace[count] = vector[nx - 3][j][k];
      counter++;
  }
  // YFACES
  counter = 0;
  //#pragma omp for collapse(2)
  for (int i = 1; i <= nxr; i++)
  for (int k = 1; k <= nzr; k++)
  {
      const int count = (i-1)*nzr+k-1;
      assert_eq(count,counter);
      ghostYleftFace[count] = vector[i][2][k];
      ghostYrghtFace[count] = vector[i][ny - 3][k];
      counter++;
  }
  // ZFACES
  counter = 0;
  //#pragma omp for collapse(2)
  for (int i = 1; i <= nxr; i++)
  for (int j = 1; j <= nyr; j++)
  {
      const int count = (i-1)*nyr+j-1;
      assert_eq(count,counter);
      ghostZleftFace[count] = vector[i][j][2];
      ghostZrghtFace[count] = vector[i][j][nz - 3];
      counter++;
  }
}

// / SPECIES for interpolation
/** prepare ghost cells on 6 faces for communication */
void makeCenterFace(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtFace, double *ghostXleftFace,
  double *ghostYrghtFace, double *ghostYleftFace,
  double *ghostZrghtFace, double *ghostZleftFace)
{
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;
  // XFACES 
  int counter = 0;
  for (int j = 1; j <= nyr; j++)
  for (int k = 1; k <= nzr; k++)
  {
      const int count = (j-1)*nzr+k-1;
      assert_eq(count,counter);
      ghostXleftFace[counter] = vector[1][j][k];
      ghostXrghtFace[counter] = vector[nx - 2][j][k];
      counter++;
  }
  // YFACES
  counter = 0;
  for (int i = 1; i <= nxr; i++)
  for (int k = 1; k <= nzr; k++)
  {
      const int count = (i-1)*nzr+k-1;
      assert_eq(count,counter);
      ghostYleftFace[counter] = vector[i][1][k];
      ghostYrghtFace[counter] = vector[i][ny - 2][k];
      counter++;
  }
  // ZFACES
  counter = 0;
  for (int i = 1; i <= nxr; i++)
  for (int j = 1; j <= nyr; j++)
  {
      const int count = (i-1)*nyr+j-1;
      assert_eq(count,counter);
      ghostZleftFace[counter] = vector[i][j][1];
      ghostZrghtFace[counter] = vector[i][j][nz - 2];
      counter++;
  }
}

// ////////////////////
// ///////////////////
// EDGES
// ///////////////////
// ///////////////////

/** prepare ghost cell Edge Z for communication: these are communicated in Y direction */
void makeNodeEdgeZ(int nx, int ny, int nz,
  double *faceXleft, double *faceXrght,
  double *ghostXrghtYrghtZsameEdge, double *ghostXleftYleftZsameEdge,
  double *ghostXrghtYleftZsameEdge, double *ghostXleftYrghtZsameEdge)
{
  const int nxr = nx - 2;
  const int nyr = ny - 2;
  const int nzr = nz - 2;
  int counter = 0;
  int counterLeft = 0;
  int counterrght = 0;
  for (int j = 1; j <= nyr; j++)
  for (int k = 1; k <= nzr; k++)
  {
      // YLEFT
      if (j == 1)
      {
        const int count = k-1;
        const int countLeft = k-1;
        assert_eq(count,counter);
        assert_eq(countLeft,counterLeft);
        ghostXleftYleftZsameEdge[counterLeft] = faceXleft[counter];
        ghostXrghtYleftZsameEdge[counterLeft] = faceXrght[counter];
        counterLeft++;
      }
      // Yrght 
      if (j == nyr)
      {
        const int count = (nyr-1)*nzr+k-1;
        const int countrght = k-1;
        assert_eq(count,counter);
        assert_eq(countrght,counterrght);
        ghostXleftYrghtZsameEdge[counterrght] = faceXleft[counter];
        ghostXrghtYrghtZsameEdge[counterrght] = faceXrght[counter];
        counterrght++;
      }
      counter++;
    }
}



// /
// 
// Y EDGE
/** prepare ghost cell Edge Y for communication: these are communicate: these are communicated in X direction */
void makeNodeEdgeY(int nx, int ny, int nz,
  double *faceZleft, double *faceZrght,
  double *ghostXrghtYsameZrghtEdge, double *ghostXleftYsameZleftEdge,
  double *ghostXleftYsameZrghtEdge, double *ghostXrghtYsameZleftEdge)
{
  const int nxr = nx - 2;
  const int nyr = ny - 2;
  const int nzr = nz - 2;
  int counter = 0;
  int counterLeft = 0;
  int counterRght = 0;
  for (int i = 1; i <= nxr; i++)
  for (int j = 1; j <= nyr; j++)
  {
      // XLEFT
      if (i == 1)
      {
        const int countLeft = j-1;
        const int count = j-1;
        assert_eq(countLeft,counterLeft);
        assert_eq(count,counter);
        ghostXleftYsameZleftEdge[countLeft] = faceZleft[count];
        ghostXleftYsameZrghtEdge[countLeft] = faceZrght[count];
        counterLeft++;
      }
      // Xrght 
      if (i == nxr)
      {
        const int countRght = j-1;
        const int count = (nxr-1)*nyr+j-1;
        assert_eq(countRght,counterRght);
        assert_eq(count,counter);
        ghostXrghtYsameZleftEdge[countRght] = faceZleft[count];
        ghostXrghtYsameZrghtEdge[countRght] = faceZrght[count];
        counterRght++;
      }
      counter++;
    }
}



// 
// 
// X EDGE
/** prepare ghost cell Edge X for communication: these are communicated in  Z direction*/
void makeNodeEdgeX(int nx, int ny, int nz,
  double *faceYleft, double *faceYrght,
  double *ghostXsameYrghtZrghtEdge, double *ghostXsameYleftZleftEdge,
  double *ghostXsameYleftZrghtEdge, double *ghostXsameYrghtZleftEdge)
{
  const int nxr = nx - 2;
  const int nyr = ny - 2;
  const int nzr = nz - 2;
  int counter = 0;
  int counterLeft = 0;
  int counterrght = 0;
  for (int i = 1; i <= nxr; i++)
  for (int k = 1; k <= nzr; k++)
  {
      // ZLEFT
      if (k == 1)
      {
        const int count = (i-1)*nzr+k-1;
        const int countLeft = i-1;
        assert_eq(count,counter);
        assert_eq(countLeft,counterLeft);
        ghostXsameYleftZleftEdge[counterLeft] = faceYleft[counter];
        ghostXsameYrghtZleftEdge[counterLeft] = faceYrght[counter];
        counterLeft++;
      }
      // Zrght 
      if (k == nzr)
      {
        const int count = (i-1)*nzr+k-1;
        const int countrght = i-1;
        assert_eq(count,counter);
        assert_eq(countrght,counterrght);
        ghostXsameYleftZrghtEdge[counterrght] = faceYleft[counter];
        ghostXsameYrghtZrghtEdge[counterrght] = faceYrght[counter];
        counterrght++;
      }
      counter++;
    }
}


// /////////////////////////////
// ////////////////////////////
// CORNER
// ///////////////////////////
// ///////////////////////////
/** prepare ghost cell Edge X for communication */
void makeNodeCorner(int nx, int ny, int nz, double *ghostXsameYrghtZrghtEdge, double *ghostXsameYleftZleftEdge, double *ghostXsameYleftZrghtEdge, double *ghostXsameYrghtZleftEdge, double *ghostXrghtYrghtZrghtCorner, double *ghostXleftYrghtZrghtCorner, double *ghostXrghtYleftZrghtCorner, double *ghostXleftYleftZrghtCorner, double *ghostXrghtYrghtZleftCorner, double *ghostXleftYrghtZleftCorner, double *ghostXrghtYleftZleftCorner, double *ghostXleftYleftZleftCorner) {
  *ghostXleftYrghtZrghtCorner = ghostXsameYrghtZrghtEdge[0];
  *ghostXrghtYrghtZrghtCorner = ghostXsameYrghtZrghtEdge[nx - 3];

  *ghostXleftYleftZrghtCorner = ghostXsameYleftZrghtEdge[0];
  *ghostXrghtYleftZrghtCorner = ghostXsameYleftZrghtEdge[nx - 3];

  *ghostXleftYrghtZleftCorner = ghostXsameYrghtZleftEdge[0];
  *ghostXrghtYrghtZleftCorner = ghostXsameYrghtZleftEdge[nx - 3];

  *ghostXleftYleftZleftCorner = ghostXsameYleftZleftEdge[0];
  *ghostXrghtYleftZleftCorner = ghostXsameYleftZleftEdge[nx - 3];
}
// /////////////////////////////
// ////////////////////////////
// PARSE
// ////////////////////////////
// ////////////////////////////
/** insert the ghost cells in the 3D physical vector */
void parseFace(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtFace, double *ghostXleftFace,
  double *ghostYrghtFace, double *ghostYleftFace,
  double *ghostZrghtFace, double *ghostZleftFace)
{
  const int nxr = nx - 2;
  const int nyr = ny - 2;
  const int nzr = nz - 2;
  // XFACES 
  int counter = 0;
  for (int j = 1; j <= nyr; j++)
  for (int k = 1; k <= nzr; k++)
  {
      const int count = (j-1)*nzr+k-1;
      assert_eq(count,counter);
      vector[0][j][k] = ghostXleftFace[counter];
      vector[nx - 1][j][k] = ghostXrghtFace[counter];
      counter++;
  }
  // YFACES
  counter = 0;
  for (int i = 1; i <= nxr; i++)
  for (int k = 1; k <= nzr; k++)
  {
      const int count = (i-1)*nzr+k-1;
      assert_eq(count,counter);
      vector[i][0][k] = ghostYleftFace[counter];
      vector[i][ny - 1][k] = ghostYrghtFace[counter];
      counter++;
  }
  // ZFACES
  counter = 0;
  for (int i = 1; i <= nxr; i++)
  for (int j = 1; j <= nyr; j++)
  {
      const int count = (i-1)*nyr+j-1;
      assert_eq(count,counter);
      vector[i][j][0] = ghostZleftFace[counter];
      vector[i][j][nz - 1] = ghostZrghtFace[counter];
      counter++;
  }
}

/** add the values of ghost cells faces to the 3D physical vector */
void addFace(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtFace, double *ghostXleftFace,
  double *ghostYrghtFace, double *ghostYleftFace,
  double *ghostZrghtFace, double *ghostZleftFace,
  VirtualTopology3D * vct)
{
  const int nxr = nx-2;
  const int nyr = ny-2;
  const int nzr = nz-2;

  int counter;
  // Xrght
  if (vct->hasXrghtNeighbor())
  {
    counter = 0;
    for (int j = 1; j <= nyr; j++)
    for (int k = 1; k <= nzr; k++)
    {
        const int count = (j-1)*nzr+k-1;
        assert_eq(count,counter);
        vector[nx - 2][j][k] += ghostXrghtFace[counter];
        counter++;
    }
  }
  // XLEFT
  if (vct->hasXleftNeighbor())
  {
    counter = 0;
    for (int j = 1; j <= nyr; j++)
    for (int k = 1; k <= nzr; k++)
    {
        const int count = (j-1)*nzr+k-1;
        assert_eq(count,counter);
        vector[1][j][k] += ghostXleftFace[counter];
        counter++;
    }
  }

  // Yrght
  if (vct->hasYrghtNeighbor())
  {
    counter = 0;
    for (int i = 1; i <= nxr; i++)
    for (int k = 1; k <= nzr; k++)
    {
        const int count = (i-1)*nzr+k-1;
        assert_eq(count,counter);
        vector[i][ny - 2][k] += ghostYrghtFace[counter];
        counter++;
    }
  }
  // Y LEFT
  if (vct->hasYleftNeighbor())
  {
    counter = 0;
    for (int i = 1; i <= nxr; i++)
    for (int k = 1; k <= nzr; k++)
    {
        const int count = (i-1)*nzr+k-1;
        assert_eq(count,counter);
        vector[i][1][k] += ghostYleftFace[counter];
        counter++;
    }
  }
  // Zrght
  if (vct->hasZrghtNeighbor())
  {
    counter = 0;
    for (int i = 1; i <= nxr; i++)
    for (int j = 1; j <= nyr; j++)
    {
        const int count = (i-1)*nyr+j-1;
        assert_eq(count,counter);
        vector[i][j][nz - 2] += ghostZrghtFace[counter];
        counter++;
    }
  }
  // ZLEFT
  if (vct->hasZleftNeighbor())
  {
    counter = 0;
    for (int i = 1; i <= nxr; i++)
    for (int j = 1; j <= nyr; j++)
    {
        const int count = (i-1)*nyr+j-1;
        assert_eq(count,counter);
        vector[i][j][1] += ghostZleftFace[counter];
        counter++;
    }
  }
}

/** insert the ghost cells Edge Z in the 3D physical vector */
void parseEdgeZ(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtYrghtZsameEdge, double *ghostXleftYleftZsameEdge,
  double *ghostXrghtYleftZsameEdge, double *ghostXleftYrghtZsameEdge)
{
  for (int i = 1; i < (nz - 1); i++)
  {
    vector[nx - 1][ny - 1][i] = ghostXrghtYrghtZsameEdge[i - 1];
    vector[0][0][i] = ghostXleftYleftZsameEdge[i - 1];
    vector[nx - 1][0][i] = ghostXrghtYleftZsameEdge[i - 1];
    vector[0][ny - 1][i] = ghostXleftYrghtZsameEdge[i - 1];
  }
}

/** insert the ghost cells Edge Z in the 3D physical vector */
void addEdgeZ(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtYrghtZsameEdge, double *ghostXleftYleftZsameEdge,
  double *ghostXrghtYleftZsameEdge, double *ghostXleftYrghtZsameEdge,
  VirtualTopology3D * vct)
{
  if (vct->hasXrghtNeighbor() && vct->hasYrghtNeighbor()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[nx - 2][ny - 2][i] += ghostXrghtYrghtZsameEdge[i - 1];
  }
  if (vct->hasXleftNeighbor() && vct->hasYleftNeighbor()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[1][1][i] += ghostXleftYleftZsameEdge[i - 1];
  }
  if (vct->hasXrghtNeighbor() && vct->hasYleftNeighbor()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[nx - 2][1][i] += ghostXrghtYleftZsameEdge[i - 1];
  }
  if (vct->hasXleftNeighbor() && vct->hasYrghtNeighbor()) {
    for (int i = 1; i < (nz - 1); i++)
      vector[1][ny - 2][i] += ghostXleftYrghtZsameEdge[i - 1];
  }
}
/** prepare ghost cell Edge Y for communication */
void parseEdgeY(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtYsameZrghtEdge, double *ghostXleftYsameZleftEdge,
  double *ghostXleftYsameZrghtEdge, double *ghostXrghtYsameZleftEdge)
{
  for (int i = 1; i < ny - 1; i++) {
    vector[nx - 1][i][nz - 1] = ghostXrghtYsameZrghtEdge[i - 1];
    vector[0][i][0] = ghostXleftYsameZleftEdge[i - 1];
    vector[0][i][nz - 1] = ghostXleftYsameZrghtEdge[i - 1];
    vector[nx - 1][i][0] = ghostXrghtYsameZleftEdge[i - 1];
  }
}
/** add the ghost cell values Edge Y to the 3D physical vector */
void addEdgeY(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtYsameZrghtEdge, double *ghostXleftYsameZleftEdge,
  double *ghostXleftYsameZrghtEdge, double *ghostXrghtYsameZleftEdge,
  VirtualTopology3D * vct)
{
  if (vct->hasXrghtNeighbor() && vct->hasZrghtNeighbor()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[nx - 2][i][nz - 2] += ghostXrghtYsameZrghtEdge[i - 1];
  }
  if (vct->hasXleftNeighbor() && vct->hasZleftNeighbor()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[1][i][1] += ghostXleftYsameZleftEdge[i - 1];
  }
  if (vct->hasXleftNeighbor() && vct->hasZrghtNeighbor()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[1][i][nz - 2] += ghostXleftYsameZrghtEdge[i - 1];
  }
  if (vct->hasXrghtNeighbor() && vct->hasZleftNeighbor()) {
    for (int i = 1; i < (ny - 1); i++)
      vector[nx - 2][i][1] += ghostXrghtYsameZleftEdge[i - 1];
  }
}
/** insert the ghost cells Edge X in the 3D physical vector */
void parseEdgeX(int nx, int ny, int nz, double ***vector,
  double *ghostXsameYrghtZrghtEdge, double *ghostXsameYleftZleftEdge,
  double *ghostXsameYleftZrghtEdge, double *ghostXsameYrghtZleftEdge)
{
  for (int i = 1; i < (nx - 1); i++) {
    vector[i][ny - 1][nz - 1] = ghostXsameYrghtZrghtEdge[i - 1];
    vector[i][0][0] = ghostXsameYleftZleftEdge[i - 1];
    vector[i][0][nz - 1] = ghostXsameYleftZrghtEdge[i - 1];
    vector[i][ny - 1][0] = ghostXsameYrghtZleftEdge[i - 1];
  }
}
/** add the ghost values Edge X to the 3D physical vector */
void addEdgeX(int nx, int ny, int nz, double ***vector,
  double *ghostXsameYrghtZrghtEdge, double *ghostXsameYleftZleftEdge,
  double *ghostXsameYleftZrghtEdge, double *ghostXsameYrghtZleftEdge,
  VirtualTopology3D * vct)
{
  if (vct->hasYrghtNeighbor() && vct->hasZrghtNeighbor()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][ny - 2][nz - 2] += ghostXsameYrghtZrghtEdge[i - 1];
  }
  if (vct->hasYleftNeighbor() && vct->hasZleftNeighbor()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][1][1] += ghostXsameYleftZleftEdge[i - 1];
  }
  if (vct->hasYleftNeighbor() && vct->hasZrghtNeighbor()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][1][nz - 2] += ghostXsameYleftZrghtEdge[i - 1];
  }
  if (vct->hasYrghtNeighbor() && vct->hasZleftNeighbor()) {
    for (int i = 1; i < (nx - 1); i++)
      vector[i][ny - 2][1] += ghostXsameYrghtZleftEdge[i - 1];
  }
}
// /////////////////
// ////////////////
// PARSE CORNERS
// ///////////////
// //////////////

/** insert the ghost cells Edge X in the 3D physical vector */
void parseCorner(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtYrghtZrghtCorner, double *ghostXleftYrghtZrghtCorner,
  double *ghostXrghtYleftZrghtCorner, double *ghostXleftYleftZrghtCorner,
  double *ghostXrghtYrghtZleftCorner, double *ghostXleftYrghtZleftCorner,
  double *ghostXrghtYleftZleftCorner, double *ghostXleftYleftZleftCorner)
{
  vector[nx - 1][ny - 1][nz - 1] = *ghostXrghtYrghtZrghtCorner;
  vector[0][ny - 1][nz - 1] = *ghostXleftYrghtZrghtCorner;
  vector[nx - 1][0][nz - 1] = *ghostXrghtYleftZrghtCorner;
  vector[0][0][nz - 1] = *ghostXleftYleftZrghtCorner;
  vector[nx - 1][ny - 1][0] = *ghostXrghtYrghtZleftCorner;
  vector[0][ny - 1][0] = *ghostXleftYrghtZleftCorner;
  vector[nx - 1][0][0] = *ghostXrghtYleftZleftCorner;
  vector[0][0][0] = *ghostXleftYleftZleftCorner;
}
/** add ghost cells values Corners in the 3D physical vector */
void addCorner(int nx, int ny, int nz, double ***vector,
  double *ghostXrghtYrghtZrghtCorner, double *ghostXleftYrghtZrghtCorner,
  double *ghostXrghtYleftZrghtCorner, double *ghostXleftYleftZrghtCorner,
  double *ghostXrghtYrghtZleftCorner, double *ghostXleftYrghtZleftCorner,
  double *ghostXrghtYleftZleftCorner, double *ghostXleftYleftZleftCorner,
  VirtualTopology3D * vct)
{
  if (vct->hasXrghtNeighbor() && vct->hasYrghtNeighbor() && vct->hasZrghtNeighbor())
    vector[nx - 2][ny - 2][nz - 2] += *ghostXrghtYrghtZrghtCorner;
  if (vct->hasXleftNeighbor() && vct->hasYrghtNeighbor() && vct->hasZrghtNeighbor())
    vector[1][ny - 2][nz - 2] += *ghostXleftYrghtZrghtCorner;
  if (vct->hasXrghtNeighbor() && vct->hasYleftNeighbor() && vct->hasZrghtNeighbor())
    vector[nx - 2][1][nz - 2] += *ghostXrghtYleftZrghtCorner;
  if (vct->hasXleftNeighbor() && vct->hasYleftNeighbor() && vct->hasZrghtNeighbor())
    vector[1][1][nz - 2] += *ghostXleftYleftZrghtCorner;
  if (vct->hasXrghtNeighbor() && vct->hasYrghtNeighbor() && vct->hasZleftNeighbor())
    vector[nx - 2][ny - 2][1] += *ghostXrghtYrghtZleftCorner;
  if (vct->hasXleftNeighbor() && vct->hasYrghtNeighbor() && vct->hasZleftNeighbor())
    vector[1][ny - 2][1] += *ghostXleftYrghtZleftCorner;
  if (vct->hasXrghtNeighbor() && vct->hasYleftNeighbor() && vct->hasZleftNeighbor())
    vector[nx - 2][1][1] += *ghostXrghtYleftZleftCorner;
  if (vct->hasXleftNeighbor() && vct->hasYleftNeighbor() && vct->hasZleftNeighbor())
    vector[1][1][1] += *ghostXleftYleftZleftCorner;

}
