
/**Add a periodic perturbation in J exp i(kx - \omega t); deltaBoB is the ratio (Delta B / B0) **/
void Particles3D::AddPerturbationJ(double deltaBoB, double kx, double ky, double Bx_mod, double By_mod, double Bz_mod, double jx_mod, double jx_phase, double jy_mod, double jy_phase, double jz_mod, double jz_phase, double B0) {

  // rescaling of amplitudes according to deltaBoB //
  double alpha;
  alpha = deltaBoB * B0 / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);
  jx_mod *= alpha;
  jy_mod *= alpha;
  jz_mod *= alpha;
  for (int i = 0; i < getNOP(); i++) {
    fetchU(i) += jx_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jx_phase);
    fetchV(i) += jy_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jy_phase);
    fetchW(i) += jz_mod / q[i] / npcel / invVOL * cos(kx * x[i] + ky * y[i] + jz_phase);
  }
}

/** Linear delta f for bi-maxwellian plasma */
double Particles3D::delta_f(double u, double v, double w, double x, double y, double kx, double ky, double omega_re, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double theta, Field * EMf) {
  const complex < double >I(0.0, 1.0);
  const double vperp = sqrt(v * v + w * w);
  const double vpar = u;
  const double kpar = kx;
  double kperp;
  if (ky == 0.0)                // because this formula is not valid for exactly parallel
    kperp = 1e-9;
  else
    kperp = ky;
  const double om_c = qom / c * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0)) / 2 / M_PI;
  const double phi = atan2(w, v);
  const double lambda = kperp * vperp / om_c;
  const complex < double >omega(omega_re, omega_i);

  const int lmax = 5;           // sum from -lmax to lmax

  double bessel_Jn_array[lmax + 2];
  double bessel_Jn_prime_array[lmax + 1];
  complex < double >a1[2 * lmax + 1], a2[2 * lmax + 1], a3[2 * lmax + 1];
  complex < double >factor, deltaf;

  // rotation of x,y
  double temp;
  temp = x;
  x = x * cos(theta) - y * sin(theta);
  y = temp * sin(theta) + y * cos(theta);


  /** for compilation issues comment this part: PUT in the math stuff */
  // calc_bessel_Jn_seq(lambda, lmax, bessel_Jn_array, bessel_Jn_prime_array);
  factor = (kpar * vperp / omega * df0_dvpar(vpar, vperp) + (1.0 - (kpar * vpar / omega)) * df0_dvperp(vpar, vperp));
  for (register int l = -lmax; l < 0; l++) {  // negative index
    a1[l + lmax] = factor / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * 0.5 * pow(-1.0, -l) * (bessel_Jn_array[-l - 1] - bessel_Jn_array[-l + 1]);
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * pow(-1.0, -l) * bessel_Jn_array[-l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * pow(-1.0, -l) * bessel_Jn_array[-l];
  }

  for (register int l = 0; l < lmax + 1; l++) { // positive index
    a1[l + lmax] = factor / lambda * bessel_Jn_array[l];
    a1[l + lmax] *= (double) l;
    a2[l + lmax] = factor * I * bessel_Jn_prime_array[l];
    a3[l + lmax] = kperp / omega * (vpar * df0_dvperp(vpar, vperp) - vperp * df0_dvpar(vpar, vperp)) / lambda * bessel_Jn_array[l];
    a3[l + lmax] *= (double) l;
    a3[l + lmax] += df0_dvpar(vpar, vperp) * bessel_Jn_array[l];
  }

  deltaf = (0.0, 0.0);
  for (register int l = -lmax; l < lmax + 1; l++) {
    deltaf += (a3[l + lmax] * Ex_mod * exp(I * Ex_phase) + a1[l + lmax] * Ey_mod * exp(I * Ey_phase) + a2[l + lmax] * Ez_mod * exp(I * Ez_phase)) / (kpar * vpar + l * om_c - omega) * exp(-I * phi * (double) l);
  }
  deltaf *= I * qom * exp(I * lambda * sin(phi)) * exp(I * (2 * M_PI * kx * x + 2 * M_PI * ky * y));

  return (real(deltaf));
}

/** apply a linear perturbation to particle distribution */
//
void Particles3D::linear_perturbation(double deltaBoB, double kx, double ky, double angle, double omega_r, double omega_i, double Ex_mod, double Ex_phase, double Ey_mod, double Ey_phase, double Ez_mod, double Ez_phase, double Bx_mod, double Bx_phase, double By_mod, double By_phase, double Bz_mod, double Bz_phase, const_arr4_double fieldForPcls) {

  double value1 = 0.0, value2 = 0.0, max_value = 0.0, min_value = 0.0, phi, n;
  // rescaling of amplitudes according to deltaBoB //

  const double alpha = deltaBoB * sqrt(EMf->getBx(1, 1, 0) * EMf->getBx(1, 1, 0) + EMf->getBy(1, 1, 0) * EMf->getBy(1, 1, 0) + EMf->getBz(1, 1, 0) * EMf->getBz(1, 1, 0)) / sqrt(Bx_mod * Bx_mod + By_mod * By_mod + Bz_mod * Bz_mod);

  Ex_mod *= alpha;
  Ey_mod *= alpha;
  Ez_mod *= alpha;
  Bx_mod *= alpha;
  By_mod *= alpha;
  Bz_mod *= alpha;



  // find the maximum value of f=1+delta_f/f0
  for (register double vpar = -2 * uth; vpar <= 2 * uth; vpar += 0.0005)
    for (register double vperp = 1e-10; vperp <= 2 * vth; vperp += 0.0005)
      for (register double X = xstart; X <= xend; X += 2 * grid->getDX())
        for (register double Y = ystart; Y <= yend; Y += 2 * grid->getDY()) {
          value1 = 1 + delta_f(vpar, vperp, 0.0, X, Y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(vpar, vperp);

          if (value1 > max_value)
            max_value = value1;


        }



  max_value *= 3.2;
  phi = 1.48409;
  n = 2.948687;                 // security factor...
  if (ns == 1) {
    max_value *= 3.0;
    phi = -1.65858;
    n = 2.917946;
  }                             // security factor...
  cout << "max-value=" << max_value << " min-value=" << min_value << endl;

  /* initialize random generator */
  srand(vct->getCartesian_rank() + 2);

  int counter = 0;
  int total_generated = 0;
  const double q_value = (qom / fabs(qom)) * ((0.19635) / npcel) * (1.0 / invVOL);
  for (int i = 1; i < grid->getNXC() - 1; i++)
  for (int j = 1; j < grid->getNYC() - 1; j++)
  {
    const double x_factor = (dx / (npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))))) + grid->getXN(i, j, 0);
    const double y_factor = (dy / npcely) + grid->getYN(i, j, 0);
    for (int ii = 0; ii < npcelx + (int) (2 * n * (cos(2 * M_PI * 0.4125 * grid->getXN(i, j, 0) + phi))); ii++)
    for (int jj = 0; jj < npcely; jj++)
    {
      const double x = (ii + .5) * x_factor;
      const double y = (jj + .5) * y_factor;
      eprintf("this was not being initialized.");
      const double z = 0.;

      // apply rejection method in velocity space
      double u,v,w;
      bool rejected = true;
      while (rejected) {
        total_generated++;
        sample_maxwellian(
          u,v,w,
          uth,vth,wth,
          u0,v0,w0);

        // test: if rand < (1+delta_f/f0)/max_value --> accepted
        if (rand() / (double) RAND_MAX <= (1 + delta_f(u, v, w, x, y, kx, ky, omega_r, omega_i, Ex_mod, Ex_phase, Ey_mod, Ey_phase, Ez_mod, Ez_phase, angle, EMf) / f0(u, sqrt(v * v + w * w))) / max_value)
          rejected = false;
      }
      counter++;
      add_new_particle(u,v,w,q_value,x,y,z,0.);
    }
  }
  assert_eq(getNOP(),counter+1);
  // if (vct->getCartesian_rank()==0)
  cout << "Rejection method: " << getNOP() / double (total_generated) * 100 << " % of particles are accepted for species " << ns << " counter=" << counter << endl;
}

double Particles3D::df0_dvpar(double vpar, double vperp) {
  double result;
  result = -2 * (vpar - u0) / uth / uth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::df0_dvperp(double vpar, double vperp) {
  double result;
  result = -2 * (vperp) / vth / vth * exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

double Particles3D::f0(double vpar, double vperp) {
  double result;
  result = exp(-(vperp * vperp / vth / vth + (vpar - u0) * (vpar - u0) / uth / uth));
  result *= 3.92e6 / pow(M_PI, 3 / 2) / vth / vth / uth;
  return (result);
}

