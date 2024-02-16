/**
# Sessile drop

A sessile drop is a drop of liquid at rest on a solid surface. In the
absence of gravity, the shape of the drop is controlled by surface
tension only. An important parameter is the "contact angle" $\theta$ between
the solid surface and the interface. In the absence of gravity, the
drop is hemispherical and it is easy to show that the relation between
the radius of the drop $R$ and its volume $V$ is (for two-dimensional
drops)
$$
V = R^2 (\theta - \sin\theta\cos\theta)
$$

To test this relation, a drop is initialised as a half-disk (i.e. the
initial contact angle is 90$^\circ$) and the contact angle is varied
between 15$^\circ$ and 165$^\circ$. The drop oscillates and eventually relaxes
to its equilibrium position. This equilibrium is exact to within
machine accuracy. The curvature along the interface is constant.

Note that shallower angles are [not accessible yet](/src/contact.h).

~~~gnuplot Equilibrium shapes for $15^\circ \leq \theta \leq 165^\circ$
set term push
set term @SVG size 640,180
set size ratio -1
unset key
unset xtics
unset ytics
unset border
set xrange [-1.6:1.6]
set yrange [0:]
plot 'out' w l, '' u (-$1):2 w l lt 1, 0 lt -1
set term pop
~~~
*/
double theta0 = 7.5*pi/180.;
#include "grid/multigrid.h"
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase-clsvof.h"
#include "./integralc.h"

/* #include "contact.h" */
/* #include "vof.h" */
/* #include "tension.h" */

//scalar f[], * interfaces = {f};
int LEVEL = 8;
coord interface_normal (Point point, scalar c);

#undef interface_normal
#define interface_normal(point, c) interface_normal (point, c)

coord interface_normal (Point point, scalar c)
{
  coord n;
  if(x < Delta){
    n.x = cos(theta0);
    n.y = sin(theta0);
  }
  else
    n = mycs (point, c);
  return n;
}

/**
To set the contact angle, we allocate a [height-function
field](/src/heights.h) and set the contact angle boundary condition on
its tangential component. */

//vector h[];
//h.t[left] = contact_angle (theta0*pi/180.);

int main()
{
  size (2 [1]);

  /**
  We use a constant viscosity. */
  
  /* const face vector muc[] = {.1,.1}; */
  /* mu = muc; */
  mu1 = 0.25;
  mu2 = 0.25;
  rho1 = 1.;
  rho2 = 1.;

  /**
  We must associate the height function field with the VOF tracer, so
  that it is used by the relevant functions (curvature calculation in
  particular). */

  //  f.height = h;
  N = 1 << LEVEL;

  /**
  We set the surface tension coefficient and run for the range of
  contact angles. */
  const scalar sigma[] = 1.;
  d.sigmaf = sigma;
  /* f.sigma = 1.; */

  //  for (theta0 = 15.*pi/180.; theta0 <= 165.*pi/180.; theta0 += 15.*pi/180.)
  run();
}

/**
The initial drop is a quarter of a circle. */

event init (t = 0)
{
  //  fraction (f, - (sq(x) + sq(y) - sq(0.5)));
  foreach()
    d[] = sqrt (sq(x) + sq(y)) - 0.5;

}

#if 0
event logfile (i++)
{
  fprintf (fout, "%g %g\n", t, normf(u.x).max);
}

event snapshots (t += 1)
{
  p.nodump = false;
  dump();
}
#endif

/**
At equilibrium (t = 10 seems sufficient), we output the interface
shape and compute the (constant) curvature. */

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

/**
We compare $R/R_0$ to the analytical expression, with $R_0=\sqrt{V/\pi}$.

~~~gnuplot
reset
set xlabel 'Contact angle (degrees)'
set ylabel 'R/R_0'
set arrow from 15,1 to 165,1 nohead dt 2
set xtics 15,15,165
plot 1./sqrt(x/180. - sin(x*pi/180.)*cos(x*pi/180.)/pi) t 'analytical', \
  'log' u 2:3 pt 7 t 'numerical'
~~~

## See also

* [Similar test with
   Gerris](http://gerris.dalembert.upmc.fr/gerris/tests/tests/sessile.html)
*/
