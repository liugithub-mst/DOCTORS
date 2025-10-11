#ifndef LEGENDRE_H_
#define LEGENDRE_H_

#include <vector>
#include <iostream>

class Quadrature;

double factorial(double x);
double fastFactorialcpu(int x);
double doubleFactorialcpu(int x);

class Legendre
{
public:
    Legendre();
    ~Legendre();

    std::vector<float> m_table;  //FixMe! make it public float vector instead of double
    struct Result {
        double value;
        double derivative;

        Result() : value(0), derivative(0) {}
        Result(double val, double deriv) : value(val), derivative(deriv) {}

    };
    double operator()(const unsigned int l, const double mu);  // Evaluate Legendre polynomial with order less than 10
    Result Pn_recursive(const unsigned int l, const double mu); // Evaluate Legendre polynomial based on recursive relation
    double Pn_recursive_noderiv(const unsigned int l, const double mu); // VS: added similar to above routine except no derivative is computed
    void Pn_roots_weights(const unsigned int l); // Find roots and weight of lth order Legendre polynomial
    inline const std::vector<double> & getWeight() const { return m_weight; }
    inline const std::vector<double> & getRoot() const { return m_root; }
    double table(const unsigned int ia1, const unsigned int ia2, const unsigned int il);
    void precompute(const Quadrature* quad, const unsigned int pn);

protected:
    bool m_precomputed;
    unsigned int m_angles;
    unsigned int m_pn;
    unsigned int m_ia1jmp;
    unsigned int m_ia2jmp;
    std::vector<double> m_root;
    std::vector<double> m_weight;

};

class AssocLegendre
{
public:
    AssocLegendre();
    ~AssocLegendre();

    float operator()(const int l, const int m, const float x);
};

class SphericalHarmonic
{
public:
    SphericalHarmonic();
    ~SphericalHarmonic();

    float normConst(const int l, const int m);
    //float operator()(const int l, const int m, const float theta, const float phi);
    float ylm_e(const int l, const int m, const float theta, const float phi);
    float ylm_o(const int l, const int m, const float theta, const float phi);
    float yl0(const int l, const float theta);

protected:
    AssocLegendre m_assoc;
};

#endif // LEGENDRE_H_
