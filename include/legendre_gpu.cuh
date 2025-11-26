
#ifndef LEGENDRE_GPU_CUH
#define LEGENDRE_GPU_CUH

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

// Functions used for sperical harmonic expansion
__device__ inline double fastFactorial(int x)
{
    // Factorial table up to 13!
    double facts[] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0,
                      5040.0, 40320.0, 362880.0, 3628800.0,
                      39916800.0, 479001600.0, 6227020800.0};
    return facts[x];
}
__device__ inline float doubleFactorial_rec(int x)
{
    if (x < 2)
        return 1.0;

    return x * doubleFactorial_rec(x - 2);
}

// Non-recursive implementation of above routine
__device__ inline double doubleFactorial(int x)
{
  double result = 1.0;
  
  if(x < 2)
  {
    return result;
  }
  
  for(int i = 2; i <= x; i++)
  {
    result *= (double)i;
  }  

  return result;

}

__device__ inline double normConst(const int l, const int m)
{
    // Costant used in Sperical Harmonic expansion, m cannot be zero
    if (m == 0)
    {
      printf("Got a m=0 case in the normConst which shouldn't handle it");
    }
    double t1 = (2.0 * l + 1.0) / (2.0 * M_PI);
    double t2 = doubleFactorial(l - m) / doubleFactorial(l + m);

    return (sqrt(t1 * t2));
}

__device__ inline double computAssocLegendre(const int l, const int m, const double x)
{
    // Compute the associated Legendre polynomial without the normalization constant
    if (m < 0 || m > l || fabs(x) > 1.0) 
    {
        printf("Invalid arguments: check that 0 <= m <= l and |x| <= 1");

    }
    double pmm = 1.0;
    if (m > 0) {
        double sqrt_1_x2 = sqrt(1.0 - x * x);
        for (int i = 1; i <= m; ++i) {
            pmm *= - (2 * i - 1) * sqrt_1_x2;
        }
    }
    if (l == m)
        return pmm;

    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) {
        return pmmp1;
    }

    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }

    return pll;
}

__device__ inline float sphericalHarmonicAddition(int il, float theta_1, float phi_1, float theta_2, float phi_2)
{
  int im;
  float sum = 0;
  float x1 = cosf(theta_1);
  float x2 = cosf(theta_2);
  float phi_diff = phi_1 - phi_2;

  for(im = 1; im <= il; im++)
    sum += 2.0 * fastFactorial(il - im)/fastFactorial(il + im) * computAssocLegendre(il, im, x1) * computAssocLegendre(il, im, x2) * cosf(im*phi_diff);

  sum += computAssocLegendre(il, 0, x1) * computAssocLegendre(il, 0, x2);
  return sum;
}

//  directly compute any l-th order legendre basis polynomial using the recursion relation: (l+1)P_{l+1}(x) = (2l+1)xP_l(x) - lP_{l-1}(x)
//     l P_l(x) = (2l-1)xP_{l-1}(x) - (l-1)P_{l-2}(x)
__device__ inline float legendre_poly(int il, float x)
{
  int n;
  float temp;
  float Pn_minus = 1; // P_0
  float Pn = (il==0) ? Pn_minus : x; //P_1

  // P_n, n>1
  for(n=2; n<=il; n++)
  {
    temp = Pn;
    Pn = ((2*n-1)*x*Pn - (n-1)*Pn_minus)/n;
    Pn_minus = temp;
  }
  return Pn;
}


#endif
