#include "Grid.hpp"
#include <fstream>

using namespace std;

namespace solution1
{
    Grid::Grid()
    {
        // Считываем X.txt
        ifstream fin(_pathX);
        fin >> _N;
        _X.resize(_N);
        for (size_t i(0); i < _N; i++)
            fin >> _X[i];
        
        fin.close();

        // Считываем elems.txt
        fin.open(_pathElems);
        fin >> _Ne;
        _elems.resize(_Ne);
        for (size_t i(0); i < _Ne; i++)
        {
            _elems[i].resize(_Ne);
            fin >> _elems[i][0] >> _elems[i][1] >> _elems[i][2] >> _elems[i][3];
        }
        fin.close();

        // Считываем mat.txt
        fin.open(_pathMat);
        fin >> _Nm;
        _mat.resize(_Nm);
        for (size_t i(0); i < _Nm; i++)
            fin >> _mat[i].first >> _mat[i].second;
        fin.close();

        // Считываем b1.txt
        fin.open(_pathB1);
        size_t N = 0;
        fin >> N;
        _b1.resize(N);
        for (size_t i = 0; i < N; i++)
            fin >> _b1[i].first >> _b1[i].second;
        fin.close();

        // Считываем b2.txt
        fin.open(_pathB2);// краевое условие II рода
        fin >> N;
        _b2.resize(N);
        for (size_t i = 0; i < N; i++)        
            fin >> _b2[i].first >> _b2[i].second;
        fin.close(); 

        //in.open("b3.txt");// краевое условие III рода
        //in >> N;
        //b3.resize(N);
        //for (int i = 0; i < N; i++)
        //    in >> b3[i];
        //in.close();
    }
}