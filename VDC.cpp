# include <cmath>
# include <cstdlib>
# include <ctime>
# include <iomanip>
# include <iostream>

using namespace std;

# include "VDC.h"

double vdc ( int i )

{
  double base_inv;
  int d;
  double r;
  double s;
  int t;
//
//  Isolate the sign.
//
  if ( i < 0 )
  {
    s = -1.0;
  }
  else
  {
    s = +1.0;
  }
//
//  Work with the magnitude of I.
//
  t = abs ( i );
//
//  Carry out the computation.
//
  base_inv = 0.5;

  r = 0.0;

  while ( t != 0 )
  {
    d = ( t % 2 );
    r = r + ( double ) ( d ) * base_inv;
    base_inv = base_inv / 2.0;
    t = ( t / 2 );
  }
//
//  Recover the sign.
//
  r = r * s;

  return r;
}


double *vdc_sequence ( int i1, int i2 )
{
  double base_inv;
  int d;
  int i;
  int i3;
  int j;
  int n;
  double *r;
  double s;
  int t;

  if ( i1 <= i2 )
  {
    i3 = +1;
  }
  else
  {
    i3 = -1;
  }

  n = abs ( i2 - i1 ) + 1;

  r = new double[n];

  j = 0;
  i = i1;

  while ( 1 )
  {
//
//  Isolate the sign.
//
    if ( i < 0 )
    {
      s = -1.0;
    }
    else
    {
      s = +1.0;
    }
//
//  Work with the magnitude of I.
//
    t = abs ( i );
//
//  Carry out the computation.
//
    base_inv = 0.5;

    r[j] = 0.0;

    while ( t != 0 )
    {
      d = ( t % 2 );
      r[j] = r[j] + ( double ) ( d ) * base_inv;
      base_inv = base_inv / 2.0;
      t = ( t / 2 );
    }
//
//  Recover the sign.
//
    r[j] = r[j] * s;

    j = j + 1;

    if ( i == i2 )
    {
      break;
    }

    i = i + i3;
  }

  return r;
}


