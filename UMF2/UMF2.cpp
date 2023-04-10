#include"Header.h"
using namespace solution;

void showmenu()
{
    cout << "1 - построение сетки и вывод координат" << endl
        << "2 - подготовка начального приближения" << endl
        << "Не целочисленное значение - выход" << endl;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    Difur object;
    object.Read();// ввод данных
    object.Priblijenie();
    object.Iteration();
    //object.Newton();
    return 0;
}