#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include "vector.h"

using namespace std;

class Func
{
public:
    int N;

    Func(const int& t_N) : N(t_N) {};

    Func() : N(0) {};

    string FileOutput()
    {
        switch (N)
        {
        case(0): return "/Quad.txt";
        case(1): return "/Rosenbrock.txt";
        case(2): return "/Function.txt";
        };
    }

    inline double func(vector<double> x)
    {
        switch (N)
        {
        case(0): return 100 * (x[1] - x[0]) * (x[1] - x[0]) + (1 - x[0]) * (1 - x[0]);
        //case(0): return -x[0]*x[0] - x[1]*x[1];
        case(1): return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
        case(2): return -2 / (1 + (x[0] - 2) * (x[0] - 2) / 9 + (x[1] - 2) * (x[1] - 2)) - 1 / (1 + (x[0] - 3) * (x[0] - 3) + (x[1] - 1) * (x[1] - 1) / 9);
        };
    }

    vector<vector<double>> H(const vector<double> x)
    {
        switch (N)
        {
        case(0): return { {202 , -200} , {-200, 200} };
        //case(0): return { {-2 , 0} , {0, -2} };
        case(1): return { {-400 * x[1] + 1200 * x[0] * x[0] + 1 , -400 * x[0]} , {-400 * x[0], 200} };
        case(2): return numerical_diff_H(x);
        };
    }

    vector<vector<double>> H1(const vector<double> x)
    {
        switch (N)
        {
        case(0): return inverse(H(x));
        case(1): return inverse(H(x));
        {
            vector<vector<double>> h1 = inverse(H(x));
            if (det(h1) > 0 && h1[0][0] > 0)
                return h1;
            else
                return { { 1, 0 }, { 0, 1 } };
        }
        case(2):
        {
            vector<vector<double>> h1 = inverse(H(x));
            if (det(h1) > 0 && h1[0][0] > 0)
                return h1;
            else
                return { { 1, 0 }, { 0, 1 } };
        }
        };
    }

    vector<double> Grad(vector<double> x)
    {
        switch (N)
        {
        case(0): return { 202 * x[0] - 200 * x[1] - 2, 200 * x[1] - 200 * x[0] };
        //case(0): return { -2*x[0], -2*x[1]  };
        case(1): return { -400 * x[1] * x[0] + 400 * x[0] * x[0] * x[0] - 2 + 2 * x[0], 200 * x[1] - 200 * x[0] * x[0] };
        case(2): return numerical_diff_Grad(x);
        };
    }

    inline double Lambda(vector<double> x, const double EPS)
    {
        switch (N)
        {
        case(0): return golden(interval(x), x, EPS);
        case(1): return golden(interval(x), x, EPS);
        case(2): return golden(interval(x), x, EPS);
        };
    }


    vector<vector<double>> inverse(vector<vector<double>> H)
    {
        double DET = det(H), c;
        H[0][1] = -H[0][1]; H[1][0] = -H[1][0];
        c = H[0][0]; H[0][0] = H[1][1]; H[1][1] = c;
        return (1 / DET) * H;
    }

    vector<vector<double>> numerical_diff_H(const vector<double> x)
    {
        vector<vector<double>> H(2, vector<double>(2));
        double h = 1e-1;
        vector<double> h1 = { h, 0 }, h2 = { 0, h };
        H[0][0] = ((func(x + 2 * h1) - func(x)) / (4 * pow(h, 2))) - ((func(x) - func(x - 2 * h1)) / (4 * pow(h, 2)));
        H[1][1] = ((func(x + 2 * h2) - func(x)) / (4 * pow(h, 2))) - ((func(x) - func(x - 2 * h2)) / (4 * pow(h, 2)));
        H[0][1] = H[1][0] = ((func(x + h1 + h2) - func(x + h1 - h2)) / (4 * pow(h, 2))) - ((func(x - h1 + h2) - func(x - h1 - h2)) / (4 * pow(h, 2)));
        return H;
    }

    vector<double> numerical_diff_Grad(const vector<double> x)
    {
        double h = 1e-1;
        vector<double> H(2), h1 = { h, 0 }, h2 = { 0, h };
        H[0] = (func(x + h1) - func(x - h1)) / (2 * h);
        H[1] = (func(x + h2) - func(x - h2)) / (2 * h);
        return H;
    }

    vector<double> interval(vector<double> x)
    {
        double xk, xk1, xk_1, x0 = 0, h, DELTA = 1e-5;
        vector<double> Sk(2);
        //Sk = (-1)*H1(x) * Grad(x);
        Sk = -(1 / norm(Grad(x))) * Grad(x);
        double f = func(x + x0 * Sk);
        if (f == func(x + (x0 + DELTA) * Sk))
            return { x0, x0 + DELTA };
        else if (f == func(x + (x0 - DELTA) * Sk))
            return { x0 - DELTA, x0 };
        else
        {
            if (f > func(x + (x0 + DELTA) * Sk))
            {
                xk = x0 + DELTA;
                h = DELTA;
            }
            else if (f > func(x + (x0 - DELTA) * Sk))
            {
                xk = x0 - DELTA;
                h = -DELTA;
            }
            else
            {
                return { x0 - DELTA , x0 + DELTA };
            }

            xk_1 = x0;

            bool exit = false;
            do
            {

                h *= 2;
                xk1 = xk + h;

                if (func(x + (xk + DELTA) * Sk) > func(x + (xk1 + DELTA) * Sk))
                {
                    xk_1 = xk;
                }
                else
                    exit = true;
            } while (!exit);

            return { xk_1, xk1 };
        }
    }

    inline double golden(vector<double> ab, vector<double> x, const double& EPS)
    {
        double SQRT5 = sqrt(5);
        double x1 = ab[0] + (3 - SQRT5) / 2 * (ab[1] - ab[0]);
        double x2 = ab[0] + (SQRT5 - 1) / 2 * (ab[1] - ab[0]);
        vector<double> Sk(2);
        //Sk = (-1) * H1(x) * Grad(x);
        Sk = -(1 / norm(Grad(x))) * Grad(x);
        double f1 = func(x + (x1 * Sk));
        double f2 = func(x + (x2 * Sk));
        double a1, b1;
        int n = 0;
        for (; abs(ab[1] - ab[0]) > EPS; n++)
        {
            a1 = ab[0], b1 = ab[1];
            if (f1 < f2)
            {
                ab[1] = x2;
                x2 = x1;
                x1 = ab[0] + (3 - SQRT5) / 2 * (ab[1] - ab[0]);
                f2 = f1;
                f1 = func(x + (x1 * Sk));
            }
            else
            {
                ab[0] = x1;
                x1 = x2;
                x2 = ab[0] + (SQRT5 - 1) / 2 * (ab[1] - ab[0]);
                f1 = f2;
                f2 = func(x + (x2 * Sk));
            }
        }
        return ab[0];
    }
};