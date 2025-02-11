#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "tension.h"

#define U 1      // Wall velocity
#define H 1      // Channel half height
#define a 0.5*H  // Drop radius
#define Re 1     // Reynolds number
#define Ca 0.9   // Capillary number
#define M 1      // Viscosity ratio


// Set walls on top and bottom boundaries
u.t[bottom] = dirichlet(-U);
u.t[top] = dirichlet(U);

int main() {

  size (2.*H);
  origin (-1., -1.);
  periodic (right);
  init_grid (64);

  rho1 = 1., rho2 = 1.;
  mu1 = rho1*U*2*a/Re;
  mu2 = mu1*M;
  f.sigma = U*mu1/Ca;
  TOLERANCE = 1e-4;
  run();
}

event init (t = 0) {
  fraction (f, sq(x) + sq(y) - sq(a));
  foreach()
    u.x[] = y;
}

event vof (i++, first);

/* void mg_print (mgstats mg) */
/* { */
/*   if (mg.i > 0 && mg.resa > 0.) */
/*     printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa, */
/* 	    mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0., */
/* 	    mg.nrelax); */
/* } */

event images (i++) {
  output_ppm (f);
}

event deformation (t += 0.1) {
  double rmax = -HUGE, rmin = HUGE ;
  foreach (reduction(max:rmax) reduction(min:rmin)) 
    if (f[] > 0 && f[] < 1) {
      coord p;
      coord n = mycs (point, f);
      double alpha = plane_alpha (f[], n);
      plane_area_center (n, alpha, &p);
      double rad  = sqrt(sq(x + Delta*p.x) + sq(y + Delta*p.y)); 
      if (rad > rmax)
	rmax = rad;
      if (rad < rmin)
	rmin = rad;
    }
  double D = (rmax - rmin)/(rmax + rmin);
  fprintf (stderr, "%g %g %g %g\n", t, rmin, rmax, D);
}

event interface (t = 3.) {
  output_facets (f, stderr);
}
