#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"
#include "DAISGram.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

int main()
{
    DAISGram dais{};
    dais.load_image("/mnt/c/Users/User/Desktop/daisgram_group_5/images/dais.bmp");
    dais.sharpen().save_image("results/prova/dais_sharp.bmp");
}