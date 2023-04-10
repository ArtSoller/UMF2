#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <stdlib.h>
#include <vector>
#include <Grid.hpp>


using namespace std;

namespace solution1
{
    class Matrix
    {
    private:
        vector<double> _elems;
    public:
        Matrix(Grid _myGrid);
        size_t GetSize()
        {
            return _elems.size();
        };

    };

    class Vector
    {

    };
}

#endif