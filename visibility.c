/*----------------------------------------------------------------------------

  Copyright (c) 2018-2020 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "iio.h"

/*----------------------------------------------------------------------------*/
#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*----------------------------------------------------------------------------*/
/* PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

/*----------------------------------------------------------------------------*/
/* Label for pixels with undefined gradient. */
#define NOTDEF -2.0

/*----------------------------------------------------------------------------*/
/* fatal error, print a message to standard error and exit
 */
static void error(char * msg)
{
  fprintf(stderr,"error: %s\n",msg);
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* memory allocation, print an error and exit if fail
 */
static void * xmalloc(size_t size)
{
  void * p;
  if( size == 0 ) error("xmalloc input: zero size");
  p = malloc(size);
  if( p == NULL ) error("out of memory");
  return p;
}

/*----------------------------------------------------------------------------*/
/* Normalized angle difference between 'a' and the symmetric of 'b'
   relative to a vertical axis.
 */
static double norm_angle(double a, double b)
{
  a -= b;
  while( a <= -M_PI ) a += 2.0*M_PI;
  while( a >   M_PI ) a -= 2.0*M_PI;

  return fabs(a) / M_PI;
}

/*----------------------------------------------------------------------------*/
static double grow_region( double * image, int X, int Y, int x0, int y0,
                           int * reg_x, int * reg_y, int * reg_n, double th )
{
  int neighbor_x[4] = {  0,  0, -1,  1 };
  int neighbor_y[4] = { -1,  1,  0,  0 };
  double sum;
  int n,k;

  /* initialize region */
  sum = image[x0+y0*X];
  image[x0+y0*X] = 1.0;   /* mark pixel as already used */
  reg_x[0] = x0;          /* add pixel to region */
  reg_y[0] = y0;
  *reg_n = 1;

  /* recursively add neighbors of the that match the criterion */
  for(n=0; n < *reg_n; n++)
    for(k=0; k<4; k++)
      {
        int x = reg_x[n] + neighbor_x[k]; /* try 4-connected neighbors */
        int y = reg_y[n] + neighbor_y[k];

        if( x > 0 && x < X-1 && y > 0 && y < Y-1 && image[x+y*X] <= th )
          {
            sum += image[x+y*X];
            image[x+y*X] = 1.0;   /* mark pixel as already used */
            reg_x[*reg_n] = x;    /* add pixel to region */
            reg_y[*reg_n] = y;
            ++(*reg_n);
          }
      }

  return sum;
}

/*----------------------------------------------------------------------------*/
/* evaluates visibility in registered time series images

   it is assumed that the output array 'visible' is allocated
 */
void visibility( double * images, double * visible, int X, int Y, int N,
                 int lambda )
{
  double rho = 0.2;
  double * angles;
  double * adiff;
  int * reg_x;
  int * reg_y;
  int reg_n;
  int * done;
  double logNT,logNFA;
  int neighbor_x[4] = {  0,  0, -1,  1 };
  int neighbor_y[4] = { -1,  1,  0,  0 };
  int i,j,k,x,y;
  double V;

  /* check input */
  if( images == NULL || visible == NULL || X <= 0 || Y <= 0 || N < 2 )
    error("visibility: invalid input");

  /* get memory and initialize data */
  angles = (double *) xmalloc( X * Y * N * sizeof(double) );
  adiff = (double *) xmalloc( X * Y * sizeof(double) );
  reg_x = (int *) xmalloc( X * Y * sizeof(int) );
  reg_y = (int *) xmalloc( X * Y * sizeof(int) );
  done = (int *) xmalloc( X * Y * sizeof(int) );
  for(k=0; k<X*Y*N; k++) visible[k] = 255.0; /* initialize as not visible */

  /* compute gradient angles */
  for(i=0; i<N; i++)
  for(x=1; x<X-1; x++)
  for(y=1; y<Y-1; y++)
    {
      double dx = images[(x+1)+ y   *X + i*X*Y] - images[(x-1)+ y   *X + i*X*Y];
      double dy = images[ x   +(y+1)*X + i*X*Y] - images[ x   +(y-1)*X + i*X*Y];
      double mod = sqrt( dx*dx + dy*dy );

      if( mod <= 0.0 ) angles[x+y*X + i*X*Y] = NOTDEF;
      else             angles[x+y*X + i*X*Y] = atan2(dy,dx);
    }

  /* compare images to find matching patterns */
  for(i=0; i<N; i++)
    for(j=i+1; j<N; j++)
      {
        /* compute gradient angle difference */
        for(x=1; x<X-1; x++)
        for(y=1; y<Y-1; y++)
          {
            double a = angles[x+y*X + i*X*Y];
            double b = angles[x+y*X + j*X*Y];
            if( a != NOTDEF && b != NOTDEF ) adiff[x+y*X] = norm_angle(a,b);
            else                             adiff[x+y*X] = 1.0;
          }

        /* grow regions as potential match candidates */
        for(x=1; x<X-1; x++)
        for(y=1; y<Y-1; y++)
          if( adiff[x+y*X] <= rho )
            {
              double s = grow_region(adiff,X,Y,x,y,reg_x,reg_y,&reg_n,rho);

              /* Number of test:
                 all pair of images ~ N^2                    ->  N^2
                 each pixel is the center of possible region ->  X*Y
                 we need to test regions of size 1 to X*Y    ->  X*Y
                 the number of polyominos of size n is about
                     0.316915 * 4.0625696^n / n              -> 0.3*4.06^n / n
               */
              logNT = 2.0 * log10(N) + 2.0 * log10(X) + 2.0 * log10(Y)
                    + log10(0.316915) + (double) reg_n * log10(4.0625696)
                    - log10( (double) reg_n );

              /* NFA = NT * s^n / n!
                   log(n!) is bounded by Stirling's approximation:
                   n! >= sqrt(2pi) * n^(n+0.5) * exp(-n)
                   then, log10(NFA) <= log10(NT) + n*log10(s)
                                                 - log10(latter expansion) */
              logNFA = logNT + (double) reg_n * log10(s)
                     - 0.5 * log10(2.0 * M_PI)
                     - ( (double) reg_n + 0.5 ) * log10( (double) reg_n )
                     + (double) reg_n * log10(exp(1.0));

              /* if a meaningful match is found,
                 mark the region as visible in both images */
              if( logNFA < 0.0 )
                for(k=0; k<reg_n; k++)
                  {
                    int xx = reg_x[k];
                    int yy = reg_y[k];
                    visible[xx+yy*X + i*X*Y] = visible[xx+yy*X + j*X*Y] = 0.0;
                  }
            }
      }

  /* grain filter: remove connected regions of less than lambda pixels */
  for(i=0; i<N; i++)
  for(V=0.0; V <= 255.0; V+=255.0) /* perform grain filter for 0 and 255 */
    {
      for(k=0; k<X*Y; k++) done[k] = FALSE;

      for(x=1; x<X-1; x++)
      for(y=1; y<Y-1; y++)
        if( ! done[x+y*X] && visible[x+y*X + i*X*Y] == V )
          {
            /* initialize region */
            done[x+y*X] = TRUE;
            reg_x[0] = x;
            reg_y[0] = y;
            reg_n = 1;

            /* recursively add neighbors of the that match the criterion */
            for(k=0; k<reg_n; k++)
            for(j=0; j<4; j++)
              {
                int xx = reg_x[k] + neighbor_x[j];
                int yy = reg_y[k] + neighbor_y[j];

                if( xx > 0 && xx < X-1 && yy > 0 && yy < Y-1 &&
                    ! done[xx+yy*X] && visible[xx+yy*X + i*X*Y] == V )
                  {
                    done[xx+yy*X] = TRUE;
                    reg_x[reg_n] = xx;
                    reg_y[reg_n] = yy;
                    ++reg_n;
                  }
              }

            /* remove connected region if smaller than lambda */
            if( reg_n < lambda )
              for(k=0; k<reg_n; k++)
                visible[ reg_x[k] + reg_y[k]*X + i*X*Y ] = (V == 255.0) ? 0.0
                                                                        : 255.0;
          }
    }

  /* correct border pixels */
  for(i=0; i<N; i++)
    {
      for(x=0;x<X;x++) visible[ x + 0   *X+i*X*Y] = visible[ x + 1   *X+i*X*Y];
      for(x=0;x<X;x++) visible[ x +(Y-1)*X+i*X*Y] = visible[ x +(Y-2)*X+i*X*Y];
      for(y=0;y<Y;y++) visible[ 0 + y   *X+i*X*Y] = visible[ 1 + y   *X+i*X*Y];
      for(y=0;y<Y;y++) visible[X-1+ y   *X+i*X*Y] = visible[X-2+ y   *X+i*X*Y];
    }

  /* free memory */
  free( (void *) angles);
  free( (void *) adiff );
  free( (void *) reg_x );
  free( (void *) reg_y );
  free( (void *) done );
}

/*----------------------------------------------------------------------------*/
/*                                    main                                    */
/*----------------------------------------------------------------------------*/
int main(int argc, char ** argv)
{
  int lambda;
  double * images;
  double * mask;
  int X,Y,N;
  int XX,YY,CC;
  int n,k;

  /* usage */
  if( argc < 4 )
    error("usage: visibility <lambda> image1 image2 [image3 ... imageN]");

  /* read input */
  lambda = atoi(argv[1]);
  N = argc - 2;
  for(n=0; n<N; n++)
    {
      double * image = iio_read_image_double_split(argv[n+2], &XX, &YY, &CC);

      if( XX < 1 || YY < 1 ) error("invalid image size");
      if( CC > 1 ) error("only single channel images handled");
      if( n>0 && ( XX!=X || YY!=Y ) ) error("images must have the same size");

      /* get memory */
      if( n == 0 )
        {
          X = XX;
          Y = YY;
          images = (double *) xmalloc( X * Y * N * sizeof(double) );
          mask   = (double *) xmalloc( X * Y * N * sizeof(double) );
        }

      /* compute copy image */
      for(k=0; k<X*Y; k++)
        images[k + n*X*Y] = image[k];

      free( (void *) image );
    }

  /* call visibility detector */
  visibility(images,mask,X,Y,N,lambda);

  /* write cloud masks */
  for(n=0; n<N; n++)
    {
      char filename[512];
      sprintf(filename,"%03d.png",n);
      iio_write_image_double_vec(filename, mask + n*X*Y, X, Y, 1);
    }

  /* free memory */
  free( (void *) images );
  free( (void *) mask );

  return EXIT_SUCCESS;
}
/*----------------------------------------------------------------------------*/
