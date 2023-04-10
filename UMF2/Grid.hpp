#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <stdint.h>

using namespace std;

namespace solution1
{
    class Grid
    {
    private:
        const string _pathX = "X.txt";
        const string _pathElems = "elems.txt";
        const string _pathMat = "mat.txt";
        const string _pathB1 = "b1.txt";
        const string _pathB2 = "b2.txt";
        const string _pathB3 = "b3.txt";
        size_t _N, _Ne, _Nm;
        vector<double> _X;
        vector<vector<uint16_t> > _elems;
        vector<pair<int, double> > _mat, _b1;
        vector<pair<int, int> > _b2;
        vector<int> _b3;
    public:
        Grid();
    };
}
#endif