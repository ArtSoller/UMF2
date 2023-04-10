#include "Header.h"

namespace solution
{
    void Difur::Read()
    {
        ifstream in("X.txt");// координаты (x)
        in >> Nn;
        q.resize(Nn);
        qn.resize(Nn);
        X.resize(Nn);
        di.resize(Nn);
        B.resize(Nn);
        gl.resize(Nn);
        gu.resize(Nn);
        for (int i = 0; i < Nn; i++)
            in >> X[i];
        in.close();
        in.open("elems.txt");
        in >> Ne;
        elems.resize(Ne);
        for (int i = 0; i < Ne; i++)
        {
            elems[i].resize(4);
            in >> elems[i][0] >> elems[i][1] >> elems[i][2] >> elems[i][3];
        }
        in.close();
        int N;
        in.open("mat.txt");
        in >> N;
        mat.resize(N);
        for (int i = 0; i < N; i++)
            in >> mat[i].first >> mat[i].second;
        in.close();
        in.open("b1.txt");// краевое условие I рода
        in >> N;
        b1.resize(N);
        for (int i = 0; i < N; i++)
            in >> b1[i].first >> b1[i].second;
        in.close();
        in.open("b2.txt");// краевое условие II рода
        in >> N;
        b2.resize(N);
        for (int i = 0; i < N; i++)        
            in >> b2[i].first >> b2[i].second;
        in.close(); 
        //in.open("b3.txt");// краевое условие III рода
        //in >> N;
        //b3.resize(N);
        //for (int i = 0; i < N; i++)
        //    in >> b3[i];
        //in.close();
    }      

    void Difur::GaMbo()
    {
        vector<vector<double>> G, M;// локальная матрица жесткости
        vector<double> b;// локальный вектор правой части
        G.resize(3);
        M.resize(3);
        for (int i = 0; i < 3; i++)
            G[i].resize(3);
        for (int i = 0; i < 3; i++)
            M[i].resize(3);
        b.resize(3);
        double h, g; // шаг
        for (int i = 0; i < Ne; i++)
        {
            h = (X[elems[i][2] - 1] - X[elems[i][0] - 1]);
            g = 1/*mat[elems[i][3] - 1].second*/;
            G[0][0] = G[2][2] = 7 * (Lambda(mat[elems[i][3] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][1] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][2] - 1])) / (3 * h);
            G[0][1] = G[1][0] = G[1][2] = G[2][1] = (- 8)* (Lambda(mat[elems[i][3] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][1] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][2] - 1])) / (3 * h);
            G[0][2] = G[2][0] = (Lambda(mat[elems[i][3] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][1] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][2] - 1])) / (3 * h);
            G[1][1] = 16 * (Lambda(mat[elems[i][3] - 1].first, q[elems[i][0] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][1] - 1]) + Lambda(mat[elems[i][3] - 1].first, q[elems[i][2] - 1])) / (3 * h);

            M[0][0] = M[2][2] = g * (h * 4 / 30);
            M[0][1] = M[1][0] = M[1][2] = M[2][1] = g * (h * 2 / 30);
            M[0][2] = M[2][0] = g * (h * (-1) / 30);
            M[1][1] = g * (h * 16 / 30);


            b[0] = h / 30. * (4 * Func(X[elems[i][0] - 1]) + 2 * Func(X[elems[i][1] - 1]) - Func(X[elems[i][2] - 1]));
            b[1] = h / 30. * (2 * Func(X[elems[i][0] - 1]) + 16 * Func(X[elems[i][1] - 1]) + 2 * Func(X[elems[i][2] - 1]));
            b[2] = h / 30. * (-Func(X[elems[i][0] - 1]) + 2 * Func(X[elems[i][1] - 1]) + 4 * Func(X[elems[i][2] - 1]));

            for (int j = 0; j < 3; j++)//Заполняем глобальный вектор и главную диагональ глобальной матрицы
            {
                B[elems[i][j] - 1] += b[j];
                di[elems[i][j] - 1] += G[j][j] + M[j][j];
            }

            gu[elems[i][0] - 1] += G[0][1] + M[0][1] + G[0][2] + M[0][2] + G[1][2] + M[1][2];
            gl[elems[i][0]] += G[1][0] + M[1][0] + G[2][0] + M[2][0] + G[2][1] + M[2][1];
        }
    }    

    //void Difur::Newton_GaMbo()
    //{
    //    vector<vector<double>> G, M;// локальная матрица жесткости
    //    vector<double> b;// локальный вектор правой части
    //    G.resize(3);
    //    M.resize(3);
    //    for (int i = 0; i < 3; i++)
    //        G[i].resize(3);
    //    for (int i = 0; i < 3; i++)
    //        M[i].resize(3);
    //    b.resize(3);
    //    double h; // шаг
    //    for (int i = 0; i < Ne; i++)
    //    {
    //        h = (X[elems[i][2] - 1] - X[elems[i][0] - 1]);
    //        l = mat[elems[i][9] - 1].first;
    //        g = mat[elems[i][9] - 1].second;

    //        G[0][0] = G[2][2] = mat[elems[i][3] - 1].first * (1 / h * 7 / 3);
    //        G[0][1] = G[1][0] = G[1][2] = G[2][1] = mat[elems[i][3] - 1].first * (1 / h * (-8) / 3);
    //        G[0][2] = G[2][0] = mat[elems[i][3] - 1].first * (1 / h * 1 / 3);
    //        G[1][1] = mat[elems[i][3] - 1].first * (1 / h * 16 / 3);

    //        M[0][0] = M[2][2] = mat[elems[i][3] - 1].second * (h * 4 / 30);
    //        M[0][1] = M[1][0] = M[1][2] = M[2][1] = mat[elems[i][3] - 1].second * (h * 2 / 30);
    //        M[0][2] = M[2][0] = mat[elems[i][3] - 1].second * (h * (-1) / 30);
    //        M[1][1] = mat[elems[i][3] - 1].second * (h * 16 / 30);


    //        for (int k = 0; k < 3; k++) //Заполняем локальный вектор
    //            for (int j = 0; j < 3; j++)
    //                b[k] += M[k][j] * Func(X[elems[i][j] - 1]);

    //        for (int j = 0; j < 3; j++)//Заполняем глобальный вектор и главную диагональ глобальной матрицы
    //        {
    //            B[elems[i][j] - 1] += b[j];
    //            di[elems[i][j] - 1] += G[j][j] + M[j][j];
    //        }

    //        gu[elems[i][0] - 1] += G[0][1] + M[0][1] + G[0][2] + M[0][2] + G[1][2] + M[1][2];
    //        gl[elems[i][0]] += G[1][0] + M[1][0] + G[2][0] + M[2][0] + G[2][1] + M[2][1];
    //    }
    //}

    void Difur::Boundary_conditions() // учет краевых условий
    {
        //for (int i = 0; i < b2.size(); i++) // краевое условие II рода
        //{
        //    B[f - 1] += h / 30 * (4 * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 2 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) - Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
        //    B[s - 1] += h / 30 * (2 * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 16 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) + 2 * Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
        //    B[t - 1] += h / 30 * ((-1) * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 2 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) + 4 * Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
        //}

        //for (int i = 0; i < b2.size(); i++) // краевое условие III рода
        //{
        //    B[f - 1] += h / 30 * (4 * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 2 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) - Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
        //    B[s - 1] += h / 30 * (2 * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 16 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) + 2 * Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
        //    B[t - 1] += h / 30 * ((-1) * Tetta(b2[i][2], XY[f - 1].first, XY[f - 1].second) + 2 * Tetta(b2[i][2], XY[s - 1].first, XY[s - 1].second) + 4 * Tetta(b2[i][2], XY[t - 1].first, XY[t - 1].second));
        //}

        for (int i = 0; i < b1.size(); i++) // краевое условие I рода
        {
            double n = b1[i].first;// Номер краевого узла
            di[n - 1] = 1;
            B[n - 1] = b1[i].second;    
            gu[n - 1] = 0;
            gl[n - 1] = 0;
        }
    }   

    void Difur::Priblijenie() // приближение
    {
        for (int i = 0; i < Nn; i++)
            q[i] = 0;
    }

    void Difur::mult_A(vector<double>& res)
    {
        for (int i = 0; i < Nn; i++)
        {
            if(i > 0)
                res[i] += gl[i] * q[i - 1];
            res[i] += di[i] * q[i];
            if(i < Nn - 1)
                res[i] += gu[i] * q[i + 1];     
        }
    }

    double Difur::norm_vector(vector<double>& vec)
    {
        double sum = 0;
        for (int i = 0; i < Nn; i++)
            sum += vec[i] * vec[i];
        return sqrt(sum);
    }

    void Difur::Factorization()
    {
        for (int i = 1; i < Nn; i++)
        {
            di[i] = di[i] - gl[i] * gu[i - 1];
            gu[i] = gu[i] / di[i];
        }
    }

    void Difur::ForwardGauss(vector<double>& res)
    {
        res[0] = B[0] / di[0];
        for (int i = 1; i < Nn; i++)
            res[i] = (B[i] - res[i - 1] * gl[i]) / di[i];
    }

    void Difur::BackwardGauss(vector<double>& res, vector<double>& f)
    {
        res[Nn - 1] = f[Nn - 1];
        for (int i = Nn - 2; i >= 0; i--)
            res[i] = (f[i] - res[i + 1] * gu[i]);
    }

    void Difur::Iteration()
    {
        vector<double> y;
        vector<double> Ax;
        vector<double> qm;
        int k = 0;//Число итераций
        y.resize(Nn);
        Ax.resize(Nn);
        qm.resize(Nn);
        GaMbo();
        Boundary_conditions();
        do
        {
            k++;
            Factorization();
            ForwardGauss(y);
            BackwardGauss(qn, y);
            for (int i = 0; i < Nn; i++)
            {
                qm[i] = qn[i] - q[i];
                q[i] = w * qn[i] + (1 - w) * q[i];
                di[i] = 0;
                B[i] = 0;
                gl[i] = 0;
                gu[i] = 0;
                Ax[i] = 0;
            }
            GaMbo();
            Boundary_conditions();
            mult_A(Ax);
            for (int i = 0; i < Nn; i++)
                Ax[i] -= B[i];
            if (k % 5 == 1)
                cout << scientific << setprecision(15) << norm_vector(Ax) / norm_vector(B) << endl;
        } while (norm_vector(Ax) / norm_vector(B) >= eps && norm_vector(qm) / norm_vector(q) >= delta && k < 1000);
        cout << "k = " << k << endl;
        for (int i = 0; i < Nn; i++)
            cout << scientific << setprecision(15) << q[i] << endl;
    }

    //void Difur::Newton()
    //{
    //    vector<double> y;
    //    vector<double> Ax;
    //    vector<double> qm;
    //    double N = 0;
    //    int k = 0;//Число итераций
    //    y.resize(Nn);
    //    Ax.resize(Nn);
    //    qm.resize(Nn);
    //    Newton_GaMbo();
    //    Boundary_conditions();
    //    do
    //    {
    //        k++;
    //        Factorization();
    //        ForwardGauss(y);
    //        BackwardGauss(qn, y);
    //        for (int i = 0; i < Nn; i++)
    //        {
    //            qm[i] = qn[i] - q[i];
    //            q[i] = w * qn[i] + (1 - w) * q[i];
    //            di[i] = 0;
    //            B[i] = 0;
    //            gl[i] = 0;
    //            gu[i] = 0;
    //            Ax[i] = 0;
    //        }
    //        Newton_GaMbo();
    //        Boundary_conditions();
    //        mult_A(Ax);
    //        for (int i = 0; i < Nn; i++)
    //            Ax[i] -= B[i];
    //        if (k % 5 == 1)
    //            cout << scientific << setprecision(15) << norm_vector(Ax) / norm_vector(B) << endl;
    //        if (abs(N - norm_vector(Ax) / norm_vector(B)) < 1e-13)
    //        {
    //            v = v / 2;
    //            for (int i = 0; i < Nn; i++)
    //            {
    //                di[i] = 0;
    //                B[i] = 0;
    //                gl[i] = 0;
    //                gu[i] = 0;
    //                Ax[i] = 0;
    //            }
    //            Newton_GaMbo();
    //            Boundary_conditions();
    //            mult_A(Ax);
    //            for (int i = 0; i < Nn; i++)
    //                Ax[i] -= B[i];
    //        }
    //        N = norm_vector(Ax) / norm_vector(B);
    //    } while (norm_vector(Ax) / norm_vector(B) >= eps && norm_vector(qm) / norm_vector(q) >= delta && k < 1000);
    //    cout << "k = " << k << endl;
    //    for (int i = 0; i < Nn; i++)
    //        cout << scientific << setprecision(15) << q[i] << endl;
    //}
}