#include "Header.hpp"
#include "Grid.hpp"
#include "Matrix.hpp"

using namespace solution;
using namespace solution1;

int main()
{
    setlocale(LC_ALL, "Russian");
    Difur object;
    object.Read();// ввод данных
    object.Priblijenie();
    object.Iteration();
    //object.Newton();

    // My code
    
    Grid myGrid;
    Matrix A(myGrid);
    
    return 0;
}