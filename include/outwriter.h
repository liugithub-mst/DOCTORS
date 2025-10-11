#ifndef OUTWRITER_H_
#define OUTWRITER_H_

#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>

#include "mesh.h"
#include "quadrature.h"
#include "xsection.h"

class Mesh;
class Quadrature;
class XSection;

class OutWriter
{
public:
    OutWriter();

    template<typename T>
    static void writeScalarFlux(std::string filename, const XSection &xs, const Mesh& mesh, const std::vector<T>& flux);

    template<typename T>
    static void writeAngularFlux(std::string filename, const XSection &xs, const Quadrature &quad, const Mesh &mesh, const std::vector<T>& flux);

    template<typename T>
    static void writeDose(std::string filename, const Mesh& mesh, const std::vector<T>& dose);

    static void writeZoneId(std::string filename, const Mesh& mesh);

    template<typename T>
    static void writeArray(std::string filename, const std::vector<T>& arry, int offset=0, int elements=-1);

    template<typename T1, typename T2>
    static void writeArray2(std::string filename, const std::vector<T1>& arry1, const std::vector<T2>& arry2);

    template<typename T1, typename T2, typename T3>
    static void writeArray3(std::string filename, const std::vector<T1>& arry1, const std::vector<T2>& arry2, const std::vector<T3>& arry3);

    static void writeFloatArrays(std::string filename, const std::vector<std::vector<float> >& arry);

};

template<typename T>
void OutWriter::writeScalarFlux(std::string filename, const XSection& xs, const Mesh& mesh, const std::vector<T>& flux)
{
    std::ofstream fout;

    fout.open(filename.c_str(), std::ios::binary);  // write output in binary file
    int groups = xs.groupCount();
    fout.write((char*) &groups, sizeof(int));
    fout.write((char*) &mesh.xElemCt, sizeof(int));
    fout.write((char*) &mesh.yElemCt, sizeof(int));
    fout.write((char*) &mesh.zElemCt, sizeof(int));

    for(unsigned int i = 0; i < mesh.xElemCt; i++)
        fout.write((char*) &mesh.xNodes[i], sizeof(float));
    for(unsigned int i = 0; i < mesh.yElemCt; i++)
        fout.write((char*) &mesh.yNodes[i], sizeof(float));
    for(unsigned int i = 0; i < mesh.zElemCt; i++)
        fout.write((char*) &mesh.zNodes[i], sizeof(float));

    if(mesh.voxelCount()*xs.groupCount() != flux.size())
        std::cout << "WARNING: OutWriter::writeScalarFlux: the mesh size did not match the data size" << std::endl;

    for(unsigned int i = 0; i < flux.size(); i++)
        fout.write((char*) &flux[i], sizeof(double));

    fout.flush();
    fout.close();
}

template<typename T>
void OutWriter::writeAngularFlux(std::string filename, const XSection &xs, const Quadrature &quad, const Mesh &mesh, const std::vector<T> &flux)
{
    std::ofstream fout;

    if(xs.groupCount() * quad.angleCount() * mesh.voxelCount() != flux.size())
        std::cout << "WARNING: OutWriter::writeScalarFlux: the mesh size did not match the data size"
                  << std::endl;

    fout.open(filename.c_str()); // Write to a text file instead of a binary file
    fout << xs.groupCount() << '\n';
    fout << quad.angleCount() << '\n';
    fout << mesh.xElemCt << '\n';
    fout << mesh.yElemCt << '\n';
    fout << mesh.zElemCt << '\n';

    for(unsigned int i = 0; i < flux.size(); i++)
        fout << flux[i] << '\n';

    /*
    fout.open(filename.c_str(), std::ios::binary);  // write output in binary file
    int groups = xs.groupCount();
    int angles = quad.angleCount();

    fout.write((char*) &groups, sizeof(int));
    fout.write((char*) &angles, sizeof(int));
    fout.write((char*) &mesh.xElemCt, sizeof(int));
    fout.write((char*) &mesh.yElemCt, sizeof(int));
    fout.write((char*) &mesh.zElemCt, sizeof(int));

    for(unsigned int i = 0; i < mesh.xElemCt; i++)
        fout.write((char*) &mesh.xNodes[i], sizeof(float));
    for(unsigned int i = 0; i < mesh.yElemCt; i++)
        fout.write((char*) &mesh.yNodes[i], sizeof(float));
    for(unsigned int i = 0; i < mesh.zElemCt; i++)
        fout.write((char*) &mesh.zNodes[i], sizeof(float));

    if(xs.groupCount() * quad.angleCount() * mesh.voxelCount() != flux.size())
        std::cout << "WARNING: OutWriter::writeScalarFlux: the mesh size did not match the data size"
                  << std::endl;
    for(unsigned int i = 0; i < flux.size(); i++)
        fout.write((char*) &flux[i], sizeof(double));
    */

    fout.flush();
    fout.close();
}

template<typename T>
void OutWriter::writeDose(std::string filename, const Mesh& mesh, const std::vector<T>& dose)
{
    std::ofstream fout;
    fout.open(filename.c_str());

    fout << 3 << '\n';
    fout << mesh.xElemCt << '\n';
    fout << mesh.yElemCt << '\n';
    fout << mesh.zElemCt << '\n';

    for(unsigned int i = 0; i < mesh.xElemCt; i++)
        fout << mesh.xNodes[i] << '\n';
    for(unsigned int i = 0; i < mesh.yElemCt; i++)
        fout << mesh.yNodes[i] << '\n';
    for(unsigned int i = 0; i < mesh.zElemCt; i++)
        fout << mesh.zNodes[i] << '\n';

    if(mesh.voxelCount() != dose.size())
        std::cout << "WARNING: OutWriter::writeDose: the mesh size did not match the data size" << std::endl;

    for(unsigned int i = 0; i < dose.size(); i++)
        fout << dose[i] << '\n';

    fout.flush();
    fout.close();
}

template<typename T>
void OutWriter::writeArray(std::string filename, const std::vector<T>& arry, int offset, int elements)
{
    std::cout << "Writing 1D data to " << filename << std::endl;
    std::ofstream fout;
    fout.open(filename.c_str());

    int stop = arry.size();
    if(elements >= 0)
        stop = std::min(stop, offset+elements);
    for(int i = offset; i < stop; i++)
        fout << arry[i] << '\n';

    fout.flush();
    fout.close();
    std::cout << "Finished writing 1D data" << std::endl;
}

template<typename T1, typename T2>
void OutWriter::writeArray2(std::string filename, const std::vector<T1>& arry1, const std::vector<T2>& arry2)
{
    std::ofstream fout;
    fout.open(filename.c_str());

    if(arry1.size() != arry2.size())
    {
        std::cout << "OutWriter::writeArray2(): Array sizes do not match!" << std::endl;

    }

    int maxtimes = MIN(arry1.size(), arry2.size());

    for(int i = 0; i < maxtimes; i++)
        fout << arry1[i] << '\t' << arry2[i] << '\n';

    fout.flush();
    fout.close();
}

template<typename T1, typename T2, typename T3>
void OutWriter::writeArray3(std::string filename, const std::vector<T1>& arry1, const std::vector<T2>& arry2, const std::vector<T3>& arry3)
{
    std::ofstream fout;
    fout.open(filename.c_str());

    if(arry1.size() != arry2.size() || arry1.size() != arry3.size())
    {
        std::cout << "OutWriter::writeArray3(): Array sizes do not match!" << std::endl;

    }

    for(int i = 0; i < arry1.size(); i++)
        fout << arry1[i] << '\t' << arry2[i] << '\t' << arry3[i] << '\n';

    fout.flush();
    fout.close();
}

#endif // OUTWRITER_H_
