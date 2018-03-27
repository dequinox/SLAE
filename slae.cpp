#include <math.h>
#include "blackbox.h"

static const double eps = 1e-9;

double alpha = 0.0;
double rsnew = 0.0;
double rsold = 0.0;
double zeta  = 0.0;
double beta  = 0.0;
unsigned int i;

int main() {
      blackbox_init();
      const unsigned int n = blackbox_size();
      double *z = (double *) malloc(n*sizeof(double));
      double *r = (double *) malloc(n*sizeof(double));
      double *b = (double *) malloc(n*sizeof(double));
      double *x = (double *) malloc(n*sizeof(double));
      double *t = (double *) malloc(n*sizeof(double));

      blackbox_rhs(b);
      for (int i = 0; i < n; ++i) {
            t[i]   = x[i] = 0.0;
            z[i]   = r[i] = b[i];
            rsold += r[i] * r[i];
      }

      int epoch = n;
      while (epoch--) {
            blackbox_mult(z, t);
            zeta = 0.0;
            for (i = 0; i < n; ++i) {
                  zeta += t[i] * z[i];
            }
            alpha = rsold / zeta;
            rsnew = 0.0;
            for (i = 0; i < n; ++i){
                  x[i]  = x[i] + alpha * z[i];
                  r[i]  = r[i] - alpha * t[i];
                  rsnew += r[i] * r[i];
            }
            if (sqrt(rsnew) <= eps) break;
            beta = rsnew / rsold;
            for (i = 0; i < n; ++i){
                  z[i] = r[i] + beta * z[i];
            }
            rsold = rsnew;
      }
      blackbox_submit(x);
      return 0;
}
