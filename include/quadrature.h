#ifndef QUADRATURE_H_
#define QUADRATURE_H_

#include <vector>
//VS: method to compute arbitrary S_n, n>=2 added.
class Quadrature
{
public:
    Quadrature();
    Quadrature(const int sn);

    std::vector<float> wt;
    std::vector<float> mu;
    std::vector<float> eta;
    std::vector<float> zi;

    inline unsigned int angleCount() const { return m_angles; }

    void loadSn(const int sn);
    void sortIntoOctants();

private:
    unsigned int m_angles;

};

#endif // QUADRATURE_H_
