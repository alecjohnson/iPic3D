
/*! Initialize the EM field with constants values or from restart */
void EMfields3D::initBEAM(double x_center, double y_center, double z_center,
  double radius)
{
  const Grid *grid = &get_grid();

  double distance;
  // initialize E and rhos on nodes
  assert (col->getRestart_status() == 0);
  // initialize
  {
    for (int i = 0; i < nxn; i++)
      for (int j = 0; j < nyn; j++)
        for (int k = 0; k < nzn; k++) {
          Ex[i][j][k] = 0.0;
          Ey[i][j][k] = 0.0;
          Ez[i][j][k] = 0.0;
          Bxn[i][j][k] = 0.0;
          Byn[i][j][k] = 0.0;
          Bzn[i][j][k] = 0.0;
          distance = (grid->getXN(i, j, k) - x_center) * (grid->getXN(i, j, k) - x_center) / (radius * radius) + (grid->getYN(i, j, k) - y_center) * (grid->getYN(i, j, k) - y_center) / (radius * radius) + (grid->getZN(i, j, k) - z_center) * (grid->getZN(i, j, k) - z_center) / (4 * radius * radius);
          // plasma
          rhons[0][i][j][k] = rhoINIT[0] / FourPI;  // initialize with constant density
          // electrons
          rhons[1][i][j][k] = rhoINIT[1] / FourPI;
          // beam
          if (distance < 1.0)
            rhons[2][i][j][k] = rhoINIT[2] / FourPI;
          else
            rhons[2][i][j][k] = 0.0;
        }
    // initialize B on centers
    for (int i = 0; i < nxc; i++)
      for (int j = 0; j < nyc; j++)
        for (int k = 0; k < nzc; k++) {
          // Magnetic field
          Bxc[i][j][k] = 0.0;
          Byc[i][j][k] = 0.0;
          Bzc[i][j][k] = 0.0;


        }
    for (int is = 0; is < ns; is++)
    {
      grid->interpN2C(rhocs, is, rhons);
      part[i].force_free(is, rhocs);
    }
  }
}


/*!Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */

/*!Add a periodic perturbation in rho exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
void EMfields3D::AddPerturbationRho(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double ne_mod, double ne_phase, double ni_mod, double ni_phase, double B0, Grid * grid) {

  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  ne_mod *= alpha;
  ni_mod *= alpha;
  // cout<<" ne="<<ne_mod<<" ni="<<ni_mod<<" alpha="<<alpha<<endl;
  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      rhons[0][i][j][0] += ne_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ne_phase);
      rhons[1][i][j][0] += ni_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + ni_phase);
    }

  for (int is = 0; is < ns; is++)
    grid->interpN2C(rhocs, is, rhons);
}

/*!Add a periodic perturbation exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) * */
void EMfields3D::AddPerturbation(double deltaBoB, double kx, double ky, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, double B0, Grid * grid) {

  double alpha;

  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;

  for (int i = 0; i < nxn; i++)
    for (int j = 0; j < nyn; j++) {
      Ex[i][j][0] += Ex_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ex_phase);
      Ey[i][j][0] += Ey_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ey_phase);
      Ez[i][j][0] += Ez_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Ez_phase);
      Bxn[i][j][0] += Bx_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bx_phase);
      Byn[i][j][0] += By_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + By_phase);
      Bzn[i][j][0] += Bz_mod * cos(kx * grid->getXN(i, j, 0) + ky * grid->getYN(i, j, 0) + Bz_phase);

    }

  // initialize B on centers
  grid->interpN2C(Bxc, Bxn);
  grid->interpN2C(Byc, Byn);
  grid->interpN2C(Bzc, Bzn);

}

