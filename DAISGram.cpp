#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

void DAISGram::load_image(string filename)
{
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for(int i = 0; i < img.get_height(); ++i)
        for(int j = 0; j < img.get_width(); ++j)
        { 
            data(i,j,0) = (float) img.red_at(j,i);
            data(i,j,1) = (float) img.green_at(j,i);    
            data(i,j,2) = (float) img.blue_at(j,i);   
        }                
}

void DAISGram::save_image(string filename)
{
    data.clamp(0,255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for(int i = 0; i < getRows(); ++i)
        for(int j = 0; j < getCols(); ++j)
            img.set_pixel(j,i,(unsigned char) data(i,j,0),(unsigned char) data(i,j,1),(unsigned char) data(i,j,2));                   

    img.write(filename);
}
 
void DAISGram::generate_random(int h, int w, int d)
{
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}

int DAISGram::getRows() {return data.rows();}
int DAISGram::getCols() {return data.cols();}
int DAISGram::getDepth() {return data.depth();}

void _swap(Tensor& T, int rhs, int lhs)
{
    Tensor tmp{T};

    for (int i = 0; i < T.rows(); ++i)
        for (int j = 0; j < T.cols(); ++j)
        {
            T(i,j,rhs) = tmp(i,j,lhs);
            T(i,j,lhs) = tmp(i,j,rhs);
        }
}

DAISGram DAISGram::warhol()
{
    if (data.depth() < 3)
        throw dimension_mismatch(); 

    Tensor RG{data}, GB{data}, RB{data};
    DAISGram result{*this};
    _swap(RG,0,1);
    _swap(GB,1,2);
    _swap(RB,0,2);

    result.data.concat(RG,1);
    BG.concat(RB,1);
    result.data.concat(BG,0);

    return result;
}

DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[])
{
    if (data.depth() != 3)
        throw dimension_mismatch();

    DAISGram result{*this};

    for (int i = 0; i < result.getRows(); ++i)
        for (int j = 0; j < result.getCols(); ++j)
        {
            int tmp = 0;

            for (int k = 0; k < result.getDepth(); ++k)
                tmp += (result.data(i,j,k) >= rgb[k] - threshold[k]) and (result.data(i,j,k) <= rgb[k] + threshold[k]);

            if (tmp == 3)
                for (int k = 0; k < result.depth(); ++k)
                    result.data(i,j,k) = bkg.data(i,j,k);
        }
            
    return result;
}