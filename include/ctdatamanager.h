#ifndef CTDATAMANAGER_H_
#define CTDATAMANAGER_H_

#include <string>

#include "mesh.h"
#include "quadrature.h"

class CtDataManager
{
public:
    CtDataManager();
    ~CtDataManager();

    Mesh* ctNumberToWater(Mesh* mesh);
    Mesh* ctNumberToHumanPhantom(Mesh* mesh);
    Mesh* ctNumberToQuickCheck(Mesh* mesh);

    Mesh* parse16(int xbins, int ybins, int zbins,
                 double xlen, double ylen, double zlen, std::string filename);

    void rotateCTvol(Mesh* meshIn, Mesh* meshOut, float rotang);

protected:
    bool m_valid;
    Mesh* m_mesh;

};
#endif // CTDATAMANAGER_H_
