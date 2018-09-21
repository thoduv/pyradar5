/*
This code comes partly from https://oceancolor.gsfc.nasa.gov/docs/ocssw/nr__spline_8c_source.html

This code is based on the cubic spline interpolation code presented in:
Numerical Recipes in C: The Art of Scientific Computing
by
William H. Press,
Brian P. Flannery,
Saul A. Teukolsky, and
William T. Vetterling .
Copyright 1988 (and 1992 for the 2nd edition)

I am assuming zero-offset arrays instead of the unit-offset arrays
suggested by the authors.  You may style me rebel or conformist
depending on your point of view.

Norman Kuring	31-Mar-1999

Norman Kuring	28-May-2015	pklo and pkhi were not getting reset in
				splint() function.  This is now fixed.
				Thx, Robert Strickland, for pointing this out.
				
Aurelien Thorette 21-Sep-2018   encapsulation and conversion to double


*/

#include <stdlib.h>

typedef struct
{
	int copy;
	int n;
	double *x;
	double *y;
	double *y2;
	int pklo;
	int pkhi;
} cspline_t;



int _spline(
            const double x[], /* x, y values of points known */
            const double y[],
            int n, /* # points in x, y  */
            double yp1, /* 1st derivative at 1, nth point, set to  */
            double ypn, /* >= 10^30 to do nat'l spline */
            double y2[] /* output, 2nd derivative to use in splint */
) {

  int	i,k;
  double	p,qn,sig,un,*u;

  u = (double*)malloc((n-1)*sizeof(double));
  if(!u) return -1;

  if(yp1 > 0.99e30)
    y2[0] = u[0] = 0.0;
  else{
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }
  for(i = 1; i < n-1; i++){
    sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig - 1.0)/p;
    u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
  }
  if(ypn > 0.99e30)
    qn = un = 0.0;
  else{
    qn = 0.5;
    un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
  }
  y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);
  for(k = n-2; k >= 0; k--){
    y2[k] = y2[k]*y2[k+1] + u[k];
  }

  free(u);
  return 0;
}

void cspline_init(cspline_t *s)
{
	s->copy=0;
	s->pklo=0;
	s->pkhi=1;	
}

int cspline_compute_nat(cspline_t *s)
{
	int i;
	
	s->y2 = malloc(sizeof(double)*s->n);
	if(!s->y2)
		return -1;
	
	i = _spline(s->x,s->y,s->n,1e30,1e30,s->y2);
	if(!i)
		free(s->y2);
	
	return i;
}

void cspline_set_xy(cspline_t *s, const double *x, const double *y, int n)
{
	s->x = (double*)x;
	s->y = (double*)y;
	s->n = n;
	s->copy = 0;
}

int cspline_alloc_xy(cspline_t *s, int n)
{
	s->x = malloc(sizeof(double)*n);
	if(!s->x)
	{
		return -1;
	}
	
	s->y = malloc(sizeof(double)*n);
	if(!s->y)
	{
		free(s->x);
		return -1;
	}
	
	s->copy = 1;
	s->n = n;
	
	return 0;
}

void cspline_free(cspline_t *s)
{
	if(!s->y2)
		free(s->y2);
	if(s->copy && s->x)
		free(s->x);
	if(s->copy && s->y)
		free(s->y);
}

/* splint applys the fit for any x: low x - high x */
int cspline_eval(cspline_t *s, double x, double *y) 
{
  const double *xa = s->x, *ya = s->y, *y2a = s->y2;
  int		klo,khi,k;
  double		h,b,a;
  int	pklo=s->pklo,pkhi=s->pkhi;

  /*
  Based on the assumption that sequential calls to this function are made
  with closely-spaced, steadily-increasing values of x, I first try using
  the same values of klo and khi as were used in the previous invocation.
  If that interval is no longer correct, I do a binary search for the
  correct interval.
  */
  if(xa[pklo] <= x && xa[pkhi] > x){
    klo = pklo;
    khi = pkhi;
  }
  else{
    klo = 0;
    khi = s->n - 1;
    while(khi - klo > 1){
      k = (khi + klo) >> 1;
      if(xa[k] > x) khi = k;
      else          klo = k;
    }
    pklo = klo;
    pkhi = khi;
  }

  h = xa[khi] - xa[klo];
  if(h == 0){
	return -1;
  }
  a = (xa[khi] - x)/h;
  b = (x - xa[klo])/h;
  *y = a*ya[klo] + b*ya[khi] +
       ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
 s->pklo=pklo;
 s->pkhi=pkhi;
  return 0;
}
