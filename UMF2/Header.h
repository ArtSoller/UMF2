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
        vector<vector<int>> elems; // ������ ������� ������� (����� ��������, ��� ��������� ����� � ��������), �������� ��������� � ������
        vector<pair<int, double>> mat, b1; // ��������(������,�����), I ������� (������ ����� ����, �������� �-���)
        vector<pair<int, int>> b2; // II ������� (������ ����� ����, ��������)
        vector<int> b3; // III ������� (������ ����� ����)
        vector<double> X, gl, gu, di, B, q, qn;//����������, ��������� ������, ������ ���� �����, ����� � ���� �����
        double w = 1.1, eps = 1e-16, delta = 1e-16, v = 1; // �������, ����� �������� � ������, �������, ������
        int Nx, Ny, Ne, Nn, max_iter = 10000; // ���-�� ��������� �� X � �� Y, ���-�� ���������, ���������� �����, ������ ���-�� ��������
    public:
        void Read(); // ��������� � �������������� ������
        void GaMbo(); // ���������� �������
        void Newton_GaMbo(); // ���������� �������
        double Func(double x) //������� ������ �����
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
        double U_betta(double q)//������� �������� �������� �������
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
        void Boundary_conditions(); // ���� ������� �������
        void Factorization();
        void ForwardGauss(vector<double>& res);
        void BackwardGauss(vector<double>& res, vector<double>& f);
        void Priblijenie();
        void Iteration();
        void Newton();
        // ���� ��������������� ������� ��� ������� ���� ����� LU+���
        void mult_A(vector<double>& res);
        double norm_vector(vector<double>& vec);
    };
}