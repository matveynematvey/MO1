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
    double maxiters = 10000, iters = 1, angle;
    string path = "newton";
    ofstream f;
    vector<double> x, xk, sk, grad;
    vector<vector<double>> H1;

    Newton()
    {
        for (int i = 0; i < 3; i++)
        {
            Func func(i);
            f.open(path + func.FileOutput());
            for (double j = 1e-3; j > 1e-7; j*=1e-1)
            {
                newton(func, j, f);
            }
            f.close();
        }
    }

    void newton(Func var, double eps, ofstream &f)
    {
        sk.resize(2); H1.resize(2, vector<double>(2)), grad.resize(2);
        iters = 1;
        x = { 10, 10 }; xk = { 10, 10 };
        f << endl << "Eps: " << eps << endl;
        f << setw(3) << "i" << setw(11) << "x" << setw(10) << "y" << setw(15) << "f(x,y)" << setw(25) << "(s1,s2)" << setw(15) << "lambda" 
            << setw(33) << "x2-x1, y2-y1, f2-f1" << setw(15) << "angle" << setw(30) << "grad, H1" << endl;
        do
        {
            H1 = var.H1(xk);                                        //Hessian matrix
            grad = var.Grad(xk);                                    //Gradient
            sk = (-1)*(H1 * grad);                                       //Descent direction
            angle = (-1) * (xk * sk) / (norm(xk) + norm((-1)*sk));  //angle between x and sk
            double lambda = var.Lambda(xk, eps);                    //length parameter
            f << setw(3) << iters << setw(11) << xk[0]  << setw(10) << xk[1] << setw(15) << var.func(xk) << setw(13) << (-1)*sk[0] << "," << setw(11) << (-1)*sk[1] << setw(15) << lambda
                << setw(10) << abs(x[0] - xk[0]) <<  setw(11) << abs(x[1] - xk[1]) << "," << setw(11) << abs(var.func(xk) - var.func(x))
                << setw(15) << angle << setw(5) << "{" << grad[0] << "," << grad[1] << "}" << ", " << setw(5) <<  "{ { " << H1[0][0] << "," << H1[0][1] << " }," << " {" << H1[1][0] << ", " << H1[1][1] << "} }" << endl;
            x = xk;
            xk = x + (lambda * sk);
        } while ((abs(var.func(xk) - var.func(x)) > eps) && ((abs(xk[0] - x[0]) > eps) || (abs(xk[1] - x[1]) > eps)) && maxiters > iters++);
    }
};