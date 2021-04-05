#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include "function.h"

using namespace std;

class Pirson
{
public:
    double len, x, x1, x2, f1, f2, iters = 1, eps = 1, a, b, delta;
    string path = "interval";
    ofstream f;

    Pirson()
    {
        Func func;
    }
};
