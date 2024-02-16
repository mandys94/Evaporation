double theta0 = 7.5*pi/180.;
//#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "./integralc.h"

int LEVEL = 9;

vector nrmm[];
coord interface_normal (Point point, scalar c);

#undef interface_normal
#define interface_normal(point, c) interface_normal (point, c)

coord interface_normal (Point point, scalar c)
{
  coord n;
  if(x < Delta){
    n.x = -cos(theta0);
    n.y = -sin(theta0);
  }
  else
    n = mycs (point, c);
  return n;
}

int main()
{
  size (2 [1]);
  
  mu1 = 0.25;
  mu2 = 0.25;
  rho1 = 1.;
  rho2 = 1.;

  N = 1 << LEVEL;

  const scalar sigma[] = 5.;
  d.sigmaf = sigma;
  run();
}

event init (t = 0)
{
  foreach()
    d[] = sqrt (sq(x - 10.*cos(pi - theta0)) + sq(y)) - 10.;

}

event logfile (i += 10)
{
  scalar unorm[];
  foreach()
    unorm[] = sqrt(sq(u.x[]) + sq(u.y[]));
  
  fprintf (ferr, "%g %g\n", t, normf(unorm).avg);
}

event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}


event end (t = 10)
{
  output_facets (f, stdout);
  
  scalar kappa[];
  foreach()
    if(f[] > 1.e-6 && f[] < 1. - 1.e-6)
      kappa[] = distance_curvature (point,d);
    else
      kappa[] = nodata;
  //  curvature (f, kappa);
  stats s = statsf (kappa);
  double R = s.volume/s.sum, V = 2.*statsf(f).sum;
  fprintf (stderr, "%d %g %.5g %.3g\n", N, theta0, R/sqrt(V/pi), s.stddev);
}

#if TREE
event adapt (i++) {
  adapt_wavelet ({f,u}, (double[]){1e-2, 1e-4, 1e-4}, LEVEL);
}
#endif

