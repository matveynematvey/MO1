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
        case(1): return 100 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) + (1 - x[0]) * (1 - x[0]);
        case(2): return 2 / (1 + (x[0] - 2) * (x[0] - 2) / 9 + (x[1] - 2) * (x[1] - 2)) + 1 / (1 + (x[0] - 3) * (x[0] - 3) + (x[1] - 1) * (x[1] - 1) / 9);
        };
    }

    vector<vector<double>> H(vector<double> x)
    {
        switch (N)
        {
        case(0): return { {202 , -200} , {-200, 200} };
        case(1): return { {-400 * x[1] + 1200 * x[0] * x[0] + 1 , -400 * x[0]} , {-400 * x[0], 200} };
        case(2): return { {202 , -200} , {-200, 200} };
        };
    }

    vector<vector<double>> H1(vector<double> x)
    {
        switch (N)
        {
        case(0): return { {1 , 1} , {1, 1.01} };
        case(1):
        {
            double C = 200 * (-400 * x[1] + 400 * x[0] * x[0] + 1);
            vector<vector<double>> h1 = { {200, 400 * x[0]} , {400 * x[0], -400 * x[1] + 1200 * x[0] * x[0] + 1 } };
            if (det(h1) > 0)
                return (1 / C) * h1;
            else
                return { { 1, 1 }, { 1, 1 } };
        }
        case(2): return { {202 , -200} , {-200, 200} };
        };
    }

    vector<double> Grad(vector<double> x)
    {
        switch (N)
        {
        case(0): return { 202 * x[0] - 200 * x[1] - 2, 200 * x[1] - 200 * x[0] };
        case(1): return { -400 * x[1] * x[0] + 400 * x[0] * x[0] * x[0] - 2 + 2 * x[0], 200 * x[1] - 200 * x[0] * x[0] };
        case(2): return { 202 * x[0] - 200 * x[1] - 2, 200 * x[1] - 200 * x[0] };
        };
    }
};


