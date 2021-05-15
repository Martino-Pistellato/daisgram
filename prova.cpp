#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.cpp"
#include "DAISGram.cpp"
#include "libbmp.cpp"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

int main()
{   
    DAISGram dais{}; 
    dais.load_image("images/dais.bmp");
    dais.emboss().save_image("results/prova/dais_equalized.bmp");

    /*DAISGram dais{};
    DAISGram bkg{}; //Codice per il greenscreen di matrix
    dais.load_image("images/greenscreen/gs_2.bmp");
    bkg.load_image("images/greenscreen/gs_2_bkg.bmp");
    int rgb[3]={80,180,80};
    float treshold[3]={80,100,80};
    dais.greenscreen(bkg,&rgb[0],&treshold[0]).save_image("results/prova/dais_green.bmp");*/

    /*DAISGram dais{}; //Codice per il greenscreen di sebastiano
    DAISGram bkg{};
    dais.load_image("images/greenscreen/gs_4.bmp");
    bkg.load_image("images/greenscreen/gs_4_bkg.bmp");
    int rgb[3]={235,235,235};
    float treshold[3]={60,60,60};
    dais.greenscreen(bkg,&rgb[0],&treshold[0]).save_image("results/prova/dais_2_green.bmp");*/
}