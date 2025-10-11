
#ifndef LEGENDRE_GPU_CUH
#define LEGENDRE_GPU_CUH

// XL: Functions used for sperical harmonic expansion
__device__ inline float fastFactorial(int x)
{
    // Factorial table up to 13!
    float facts[] = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0,
                      5040.0, 40320.0, 362880.0, 3628800.0,
                      39916800.0, 479001600.0, 6227020800.0};
    return facts[x];
}

// VS: Since this involves recursion, it gives warnings during compilation associated with recursion:
//     "ptxas warning : Stack size for entry function <fn> ..." where fn is any function that invokes the below recursive __device__ function
//     while these warnings are pretty "normal", a non-recurive implemnattion will supress them and perhaps be more efficient.
__device__ inline float doubleFactorial_rec(int x)
{
    if (x < 2)
        return 1.0;

    return x * doubleFactorial_rec(x - 2);
}

// VS: Non-recursive implementation of above routine
__device__ inline float doubleFactorial(int x)
{
  if(x<2)
    return 1;

  int prod=1;
  int x_in = 2-(x%2);

  while(x_in <= x)
  {
    prod *= x_in;
    x_in += 2;
  }
  return prod;
}

__device__ inline float normConst(const int l, const int m)
{
    // Costant used in Sperical Harmonic expansion, m cannot be zero
    float t1 = (2.0 * l + 1.0) / (2.0 * M_PI);
    float t2 = doubleFactorial(l - m) / doubleFactorial(l + m);

    return (sqrt(t1 * t2));
}

__device__ inline float computAssocLegendre(const int l, const int m, const float x)
{
    // Compute the associated Legendre polynomial without the normalization constant
    if (m < 0 || m > l || std::abs(x) > 1.0) 
    {
        printf("Invalid arguments: check that 0 <= m <= l and |x| <= 1");

    }
    float pmm = 1.0;
    if (m > 0) {
        float sqrt_1_x2 = std::sqrt(1.0 - x * x);
        for (int i = 1; i <= m; ++i) {
            pmm *= - (2 * i - 1) * sqrt_1_x2;
        }
    }
    if (l == m)
        return pmm;

    float pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) {
        return pmmp1;
    }

    float pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll) {
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm = pmmp1;
        pmmp1 = pll;
    }

    return pll;
}


// VS: The below computes any l-th legendre basis polynomial P_l(\omega_1 . \omega_2) where \omega_i, i=1,2, is a unit vector, using associated legendre polynomials P_l^m, m=1,...,l.
// VS: We use spherical coordinate representation for direction:
//     Unit vector \omega = cos(\theta) \hat{x} + sin(\theta)cos(\phi) \hat{y} + sin(\theta)sin(\phi) \hat{z}
// We can easily show that Dot prudct between two unit vectors \omega_1 . \omega_2 = cos(\theta_1)cos(\theta_2) + sin(\theta_1)sin(\theta_2) cos(\phi_1-\phi_2)
// VS: To specially compute legendre polynomial value for an input of the form \omega_1 . \omega_2, ...
//     We use spherical harmonic addition theorem. See https://mathworld.wolfram.com/SphericalHarmonicAdditionTheorem.html
//     The below computation implements/uses the above theorem.
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

// VS: directly compute any l-th order legendre basis polynomial using the recursion relation: (l+1)P_{l+1}(x) = (2l+1)xP_l(x) - lP_{l-1}(x)
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
