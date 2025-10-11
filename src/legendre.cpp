#define _USE_MATH_DEFINES
#include <cmath>

#include <limits>
#include <iostream>

#include "legendre.h"
#include "quadrature.h"

double factorial(double x)
{
    if (x <= 1.001)
        return 1.0;

    return x * factorial(x - 1.0);
}

double fastFactorialcpu(int x)
{
    const std::vector<double> facts = { 1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0,
                                        5040.0, 40320.0, 362880.0, 3628800.0,
                                      39916800.0, 479001600.0, 6227020800.0 };
    return facts[x];
}

double doubleFactorialcpu(int x)
{
    if (x < 2)
        return 1.0;

    return x * doubleFactorialcpu(x - 2);
}

Legendre::Legendre() : m_table(), m_precomputed(false), m_angles(0), m_pn(0) {} 

Legendre::~Legendre() {}

double Legendre::operator()(const unsigned int l, const double mu)
{
    double result;
    switch (l)
    {
    case 0:
        result = 1.0;
        break;
    case 1:
        result = mu;
        break;
    case 2:
        result = 0.5 * (3 * mu * mu - 1); 
        break;
    case 3:
        result = 0.5 * (5 * pow(mu, 3.0) - 3 * mu);
        break;
    case 4:
        result = 0.125 * (35 * pow(mu, 4.0) - 30 * pow(mu, 2.0) + 3);
        break;
    case 5:
        result = 0.125 * (63 * pow(mu, 5.0) - 70 * pow(mu, 3.0) + 15 * mu);
        break;
    case 6:
        result = 0.0625 * (231 * pow(mu, 6.0) - 315 * pow(mu, 4.0) + 105 * pow(mu, 2.0) - 5);
        break;
    case 7:
        result = 0.0625 * (429 * pow(mu, 7.0) - 693 * pow(mu, 5.0) + 315 * pow(mu, 3.0) - 35 * mu);
        break;
    case 8:
        result = 0.0078125 * (6435 * pow(mu, 8.0) - 12012 * pow(mu, 6.0) + 6930 * pow(mu, 4.0) - 1260 * pow(mu, 2.0) + 35);
        break;
    case 9:
        result = 0.0078125 * (12155 * pow(mu, 9.0) - 25740 * pow(mu, 7.0) + 18018 * pow(mu, 5.0) - 4620 * pow(mu, 3.0) + 315 * mu);
        break;
    case 10:
        result = 0.00390625 * (46189 * pow(mu, 10.0) - 109395 * pow(mu, 8.0) + 90090 * pow(mu, 6.0) - 30030 * pow(mu, 4.0) + 3465 * pow(mu, 2.0) - 63);
        break;
    default:{
        std::cout << "Legendre polynomials order is greater than 10, using recursive function!" << l << std::endl;
        result = Pn_recursive_noderiv(l, mu);
        break;
        }

    };

    return static_cast<double>(result);
}

Legendre::Result Legendre::Pn_recursive(const unsigned int l, const double mu)
{
    Result result(mu, 1);  // P1 = mu;

    if(l==0)
    {
      result.value = 1;    //P0=1
      result.derivative = 0;
    }
    else
    {
      double Pn_minus_1 = 1; // P0 = 1;
      const double f = (mu != 1) ?  (1 / (mu * mu - 1)) : 0 ;
      for (int lorder = 2; lorder <= (int)l; lorder++)
      {
          const double value = ((2 * lorder - 1) * mu * result.value - (lorder - 1) * Pn_minus_1) / lorder; // Pn recursive
          result.derivative = lorder * f * (mu * value - result.value); // Pn derivative

          Pn_minus_1 = result.value;
          result.value = value;
      }
    }

    return result;
}

double Legendre::Pn_recursive_noderiv(const unsigned int l, const double mu)
{
    double Pn_minus_1 = 1; // P0 = 1;
    double value = (l==0) ? Pn_minus_1 : mu;
    double value_copy;

    for (int lorder = 2; lorder <= (int)l; lorder++)
    {
        value_copy = value;
        value = ((2 * lorder - 1) * mu * value - (lorder - 1) * Pn_minus_1) / lorder; // Pn recursive
        Pn_minus_1 = value_copy;
    }
    return value;
}


void Legendre::Pn_roots_weights(const unsigned int l)
{
    // Initialize
    const double EPSILON = 1e-15;
    m_root.resize(l);
    m_weight.resize(l);

    for (int lorder = 1; lorder <= (int)l; lorder++)
    {
        double root = cos(M_PI * (lorder - 0.25) / (l + 0.5)); // initial guess
        Result result = Pn_recursive(l, root);

        double newtonRaphsonRatio = 0.0;
        do {
            newtonRaphsonRatio = result.value / result.derivative;
            root -= newtonRaphsonRatio;
            result = Pn_recursive(l, root);
        } while (fabs(newtonRaphsonRatio) > EPSILON);

        m_root[lorder - 1] = root;
        m_weight[lorder - 1] = 2.0 / ((1 - root * root) * result.derivative * result.derivative);

    }

}

double Legendre::table(const unsigned int ia1, const unsigned int ia2, const unsigned int il)
{
    if (!m_precomputed)
    {
        std::cout << "Attempted to access a data table before it was computed!" << std::endl;
        return static_cast<double>(1.0);
    }

    if (ia1 >= m_angles || ia2 >= m_angles)
    {
        std::cout << "Attempted to access an angle that doesn't exist!  ia1="
            << ia1 << "   ia2=" << ia2 << std::endl;
        return static_cast<double>(1.0);
    }

    if (il > m_pn)
    {
        std::cout << "Attempted to access a Legendre coeff that doesn't exist!  il=" << il << std::endl;
        return static_cast<double>(1.0);
    }

    return m_table[ia1 * m_ia1jmp + ia2 * m_ia2jmp + il];
}

void Legendre::precompute(const Quadrature* quad, const unsigned int pn)
{
    if (m_precomputed)
    {
        std::cout << "Computing a table after it has already been computed!" << std::endl;
        m_table.clear();
    }

    m_angles = quad->angleCount();
    m_pn = pn;
    m_table.resize(m_angles * m_angles * (m_pn + 1));

    // These are stored for navigating the table quickly later
    m_ia1jmp = m_angles * (m_pn + 1);
    m_ia2jmp = m_pn + 1;

    for (unsigned int ia1 = 0; ia1 < m_angles; ia1++)
        for (unsigned int ia2 = 0; ia2 < m_angles; ia2++)
            for (unsigned int il = 0; il <= m_pn; il++)
            {
                m_table[ia1 * m_ia1jmp + ia2 * m_ia2jmp + il] = (*this)(il, quad->mu[ia1] * quad->mu[ia2] + quad->eta[ia1] * quad->eta[ia2] + quad->zi[ia1] * quad->zi[ia2]);
            }

    m_precomputed = true;
}

AssocLegendre::AssocLegendre() {}

AssocLegendre::~AssocLegendre() {}


float AssocLegendre::operator()(const int l, const int m, const float x)
{
    // Compute the associated Legendre polynomial without the normalization constant
    if (m < 0 || m > l || std::abs(x) > 1.0) 
    {
        std::cout << "Invalid arguments: check that 0 <= m <= l and |x| <= 1" << std::endl;

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


SphericalHarmonic::SphericalHarmonic() {}

SphericalHarmonic::~SphericalHarmonic() {}

float SphericalHarmonic::normConst(const int l, const int m)
{
    if (m == 0)
    {
        std::cout << "Got a m=0 case in the normConst which shouldn't handle it" << std::endl;
    }

    float t1 = static_cast<float>((2 * l + 1) / (2.0 * M_PI));

    float t2 = factorial(l - fabs(m)) / factorial(l + fabs(m));
    
    return (sqrt(t1 * t2));  // Normalization constant for spherical harmonics

}

float SphericalHarmonic::ylm_e(const int l, const int m, const float theta, const float phi)
{
    if (m <= 0)
    {
        std::cout << "Got a neg. m value!" << std::endl;
    }
    return normConst(l, m) * m_assoc(l, m, cos(theta)) * cos(m * phi); // Spherical harmonic defination cosine term

}

float SphericalHarmonic::ylm_o(const int l, const int m, const float theta, const float phi)
{
    if (m <= 0)
    {
        std::cout << "Got a neg m value!" << std::endl;
    }
    return normConst(l, m) * m_assoc(l, m, cos(theta)) * sin(m * phi); // Spherical harmonic defination sine term

}

float SphericalHarmonic::yl0(const int l, const float theta)
{
    // Spherical harmonic defination when m=0
    return sqrt(2*l + 1) * m_assoc(l, 0, cos(theta));  // Note:  the quandrature set is already normalized by 4Pi

}
