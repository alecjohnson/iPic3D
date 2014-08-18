
#include <mpi.h>
#include "GMRES.h"
#include "debug.h"

void GMRES(FIELD_IMAGE FunctionImage, double *xkrylov, int xkrylovlen, const double *b, int m, int max_iter, double tol, Grid * grid, VirtualTopology3D * vct, Field * field) {
  if (m > xkrylovlen) {
    if (vct->getCartesian_rank() == 0)
      cerr << "In GMRES the dimension of Krylov space(m) can't be > (length of krylov vector)/(# processors)" << endl;
    return;
  }
  bool GMRESVERBOSE = false;
  double initial_error, rho_tol, av, mu, htmp, tmp, delta = 0.001;
  double *r = new double[xkrylovlen];
  double *im = new double[xkrylovlen];
  double *v;
  double *w;


  double *s = new double[m + 1];
  double *cs = new double[m + 1];
  double *sn = new double[m + 1];
  double *y = new double[m + 1];
  int k;
  eqValue(0.0, s, m + 1);
  eqValue(0.0, cs, m + 1);
  eqValue(0.0, sn, m + 1);
  eqValue(0.0, y, m + 1);


  // allocate H for storing the results from decomposition
  double **H = newArr2(double, m + 1, m);
  for (int ii = 0; ii < m + 1; ii++)
    for (int jj = 0; jj < m; jj++)
      H[ii][jj] = 0;
  // allocate V
  double **V = newArr2(double, m+1, xkrylovlen);
  for (int ii = 0; ii < m+1; ii++)
    for (int jj = 0; jj < xkrylovlen; jj++)
      V[ii][jj] = 0;



  if (GMRESVERBOSE && vct->getCartesian_rank() == 0) {
    cout << "------------------------------------" << endl;
    cout << "-             GMRES                -" << endl;
    cout << "------------------------------------" << endl;
    cout << endl;
  }

  double normb = normP(b, xkrylovlen);
  if (normb == 0.0)
    normb = 1.0;

  for (register int itr = 0; itr < max_iter; itr++) {

    // r = b - A*x
    (field->*FunctionImage) (im, xkrylov, grid, vct);
    sub(r, b, im, xkrylovlen);
    initial_error = normP(r, xkrylovlen);

    if (itr == 0) {
      if (vct->getCartesian_rank() == 0)
        cout << "Initial residual: " << initial_error << " norm b vector (source) = " << normb << endl;
      rho_tol = initial_error * tol;

      if ((initial_error / normb) <= tol) {
        if (vct->getCartesian_rank() == 0)
          cout << "GMRES converged without iterations: initial error < tolerance" << endl;
        delete[]r;
        delete[]im;
        delete[]s;
        delete[]cs;
        delete[]sn;
        delete[]y;
        delArr2(H, m + 1);
        delArr2(V, m+1);

        return;
      }
    }

    v = V[0];
    scale(v, r, (1.0 / initial_error), xkrylovlen);
    eqValue(0.0, s, m + 1);
    s[0] = initial_error;
    k = 0;
    while (rho_tol < initial_error && k < m) {

      // w= A*V(:,k)
      v = V[k];
      w = V[k+1];
      (field->*FunctionImage) (w, v, grid, vct);
      av = normP(w, xkrylovlen);

      for (register int j = 0; j <= k; j++) {
        v = V[j];
        H[j][k] = dotP(w, v, xkrylovlen);
        addscale(-H[j][k], w, v, xkrylovlen);

      }
      H[k + 1][k] = normP(w, xkrylovlen);

      if (av + delta * H[k + 1][k] == av) {

        for (register int j = 0; j <= k; j++) {
          v = V[j];
          htmp = dotP(w, v, xkrylovlen);
          H[j][k] = H[j][k] + htmp;
          addscale(-htmp, w, v, xkrylovlen);
        }
        H[k + 1][k] = normP(w, xkrylovlen);
      }
      scale(w, (1.0 / H[k + 1][k]), xkrylovlen);

      if (0 < k) {

        for (register int j = 0; j < k; j++)
          ApplyPlaneRotation(H[j + 1][k], H[j][k], cs[j], sn[j]);

        getColumn(y, H, k, m + 1);
      }

      mu = sqrt(H[k][k] * H[k][k] + H[k + 1][k] * H[k + 1][k]);
      cs[k] = H[k][k] / mu;
      sn[k] = -H[k + 1][k] / mu;
      H[k][k] = cs[k] * H[k][k] - sn[k] * H[k + 1][k];
      H[k + 1][k] = 0.0;

      ApplyPlaneRotation(s[k + 1], s[k], cs[k], sn[k]);
      initial_error = fabs(s[k]);
      k++;
    }

    k--;
    y[k] = s[k] / H[k][k];

    for (register int i = k - 1; i >= 0; i--) {
      tmp = 0.0;
      for (register int l = i + 1; l <= k; l++)
        tmp += H[i][l] * y[l];
      y[i] = (s[i] - tmp) / H[i][i];

    }


    for (register int j = 0; j < k; j++)
    {
      const double yj = y[j];
      double* Vj = V[j];
      for (int i = 0; i < xkrylovlen; i++)
        xkrylov[i] += yj * Vj[i];
    }

    if (initial_error <= rho_tol) {
      if (vct->getCartesian_rank() == 0)
        cout << "GMRES converged at restart # " << itr << "; iteration #" << k << " with error: " << initial_error / rho_tol * tol << endl;
      delete[]r;
      delete[]im;
      delete[]s;
      delete[]cs;
      delete[]sn;
      delete[]y;
      delArr2(H, m + 1);
      delArr2(V, m + 1);
      return;
    }
    if (vct->getCartesian_rank() == 0 && GMRESVERBOSE)
      cout << "Restart: " << itr << " error: " << initial_error / rho_tol * tol << endl;

  }


  if (vct->getCartesian_rank() == 0)
    cout << "GMRES not converged !! Final error: " << initial_error / rho_tol * tol << endl;

  delete[]r;
  delete[]im;
  delete[]s;
  delete[]cs;
  delete[]sn;
  delete[]y;
  delArr2(H, m + 1);
  delArr2(V, m + 1);
  return;
}


void ApplyPlaneRotation(double &dx, double &dy, double &cs, double &sn) {
  double temp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}
