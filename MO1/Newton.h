#pragma once
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include "function.h"

using namespace std;

class Newton
{
public:
    double maxiters = 10000, iters = 1, lambda = 0.99999;
    string path = "newton";
    ofstream f;
    vector<double> x, xk, Sk;

    Newton()
    {
        for (int i = 0; i < 2; i++)
        {
            for (double j = 1e-3; j > 1e-7; j*=1e-1)
            {
                Func func(i);
                newton(func, j);
            }
        }
    }

    void newton(Func var, double eps)
    {
        Sk.resize(2);
        x = { -10, -11 }; xk = { -10, -11 };
        f.open(path + var.FileOutput());
        f << setw(3) << "i" << setw(20) << "x,y" << setw(15) << "f(x,y)" << setw(10) << "(s1,s2)" << setw(10) << "lambda" 
            << setw(15) << "x2-x1, y2-y1, f2-f1" << setw(10) << "angle" << setw(15) << "grad, H1" << endl;
        do
        {
            Sk = (-1) * var.H1(x) * var.Grad(x);
            f << setw(3) << iters << setw(10) << xk[0] << "," << setw(9) << xk[1] << setw(15) << var.func(xk) << endl;
            x = xk;
            xk = x - pow(lambda, iters) * (var.H1(x) * var.Grad(x));
        } while (abs(var.func(xk) - var.func(x) > eps || ((abs(xk[0] - x[0]) > eps) || (abs(xk[1] - x[1]) > eps))) && maxiters > iters++);
        cout << xk[0] << " " << xk[1] << " " << iters << " Eps = " << var.func(xk) << " " << var.func(x);
    }
};