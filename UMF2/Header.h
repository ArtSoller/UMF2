#pragma once
#include<iostream>
#include<fstream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

namespace solution
{
    class Difur
    {
    private:
        vector<vector<int>> elems; // второе краевое условие (номер элемента, его локальная грань и значение), элементы разбиения с узлами
        vector<pair<int, double>> mat, b1; // материал(лямбда,гамма), I краевое (глобал номер узла, значение ф-ции)
        vector<pair<int, int>> b2; // II краевое (глобал номер узла, значение)
        vector<int> b3; // III краевое (глобал номер узла)
        vector<double> X, gl, gu, di, B, q, qn;//координаты, ленточный формат, вектор прав части, текущ и след прибл
        double w = 1.1, eps = 1e-16, delta = 1e-16, v = 1; // границы, коэфф разрядки и релакс, эпсилон, дельта
        int Nx, Ny, Ne, Nn, max_iter = 10000; // кол-во разбиений по X и по Y, кол-во элементов, количество узлов, максим кол-во итераций
    public:
        void Read(); // считываем и подготавливаем данные
        void GaMbo(); // построение матрицы
        void Newton_GaMbo(); // построение матрицы
        double Func(double x) //функция правой части
        {
            return 3;
        }
        double Tetta(int n, double q)
        {
            return q + 2;
        }
        double Betta(double q)
        {
            return q;
        }
        double U_betta(double q)//функция третьего краевого условия
        {
            return q;
        }
        double dL_dq(double q)
        {
            return 1;
        }
        double Lambda(int n, double q)
        {
            return 4;
        }
        void Boundary_conditions(); // учет краевых условий
        void Factorization();
        void ForwardGauss(vector<double>& res);
        void BackwardGauss(vector<double>& res, vector<double>& f);
        void Priblijenie();
        void Iteration();
        void Newton();
        // Блок вспомагательных функций для решения СЛАУ через LU+ЛОС
        void mult_A(vector<double>& res);
        double norm_vector(vector<double>& vec);
    };
}