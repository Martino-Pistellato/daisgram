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
    /*DAISGram dais{};
    DAISGram bkg{};
    dais.load_image("images/greenscreen/gs_1.bmp");
    bkg.load_image("images/greenscreen/gs_2_bkg.bmp");
    int rgb[3]={0,128,0};
    float treshold[3]={10,20,10};
    dais.greenscreen(bkg,&rgb[0],&treshold[0]).save_image("results/prova/dais_green.bmp");*/

    DAISGram emb{};
    emb.load_image("images/flower_hires.bmp");
    emb.emboss().save_image("results/prova/dais_emboss.bmp");
}