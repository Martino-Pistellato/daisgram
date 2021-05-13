#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.cpp"
#include "DAISGram.cpp"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

int main()
{   
    DAISGram dais{};
    dais.load_image("images/dais.bmp");
    dais.warhol().save_image("results/prova/dais_warhol.bmp");
}