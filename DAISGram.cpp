#include <iostream>
#include <string>

#include "dais_exc.h"
#include "tensor.h"
#include "libbmp.h"
#include "DAISGram.h"

using namespace std;

DAISGram::DAISGram()
{
    data = Tensor();
}

DAISGram::~DAISGram()
{
    data.~Tensor();
}

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

int DAISGram::getRows() {return data.rows();}

int DAISGram::getCols() {return data.cols();}

int DAISGram::getDepth() {return data.depth();}

DAISGram DAISGram::brighten(float bright)
{
    DAISGram result{*this};

    result.data = result.data + bright;
    result.data.clamp(0,255);

    return result;
}

DAISGram DAISGram::grayscale() 
{
    DAISGram result{*this};
    float avg = 0.;

    for(int i = 0; i < getRows(); ++i)
        for(int j = 0; j < getCols(); ++j)
        {
            for(int k = 0; k < getDepth(); ++k) avg += (*this).data(i,j,k);
            for(int k = 0; k < getDepth(); ++k) result.data(i,j,k) = avg / (float)getDepth();
            avg = 0;
        }

    return result;
}

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

    result.data = result.data.concat(RG,1);
    GB = GB.concat(RB,1);
    result.data = result.data.concat(GB,0);

    return result;
}

DAISGram DAISGram::sharpen()
{
    DAISGram T{*this};
    Tensor filtro{3,3,1};

    filtro(0,0,0) = 0;
    filtro(0,1,0) = -1;
    filtro(0,2,0) = 0;
    filtro(1,0,0) = -1;
    filtro(1,1,0) = 5;
    filtro(1,2,0) = -1;
    filtro(2,0,0) = 0;
    filtro(2,1,0) = -1;
    filtro(2,2,0) = 0;

    T.data = T.data.convolve(filtro);
    T.data.clamp(0,255); 
    
    return T;
}

DAISGram DAISGram::emboss()
{
    DAISGram T{*this};
    Tensor filtro{3,3,1};

    filtro(0,0,0) = -2;
    filtro(0,1,0) = -1;
    filtro(0,2,0) = 0;
    filtro(1,0,0) = -1;
    filtro(1,1,0) = 1;
    filtro(1,2,0) = 1;
    filtro(2,0,0) = 0;
    filtro(2,1,0) = 1;
    filtro(2,2,0) = 2;

    T.data = T.data.convolve(filtro);
    T.data.clamp(0,255);
    
    return T;
}

DAISGram DAISGram::smooth(int h)    
{   //TODO risolvere con interi maggiori, tipo 7
    DAISGram T{*this};
    float c = 1 / (float)(h * h);
    Tensor filtro{3,3,1,c};
    
    T.data = T.data.convolve(filtro);
    T.data.rescale(255);
    
    return T;
}

DAISGram DAISGram::edge()
{
    DAISGram T{*this};
    Tensor filtro{3,3,1};

    filtro(0,0,0) = -1;
    filtro(0,1,0) = -1;
    filtro(0,2,0) = -1;
    filtro(1,0,0) = -1;
    filtro(1,1,0) = 8;
    filtro(1,2,0) = -1;
    filtro(2,0,0) = -1;
    filtro(2,1,0) = -1;
    filtro(2,2,0) = -1;

    T = T.grayscale();
    T.data = T.data.convolve(filtro);
    T.data.clamp(0,255);
    
    return T;
}

DAISGram DAISGram::blend(const DAISGram & rhs, float alpha)
{
    if(getCols() != rhs.data.cols() or getRows() != rhs.data.rows() or getDepth() != rhs.data.depth()) throw dimension_mismatch();
    DAISGram result{*this};

    result.data = data * alpha + rhs.data * (1 - alpha);

    return result;
}

DAISGram DAISGram::greenscreen(DAISGram & bkg, int rgb[], float threshold[]) //threshold = margine errore
{
    if (data.depth() != 3)
        throw dimension_mismatch();

    DAISGram result{*this};

    for (int i = 0; i < result.getRows(); ++i)
        for (int j = 0; j < result.getCols(); ++j)
        {
            int tmp = 0;

            for (int k = 0; k < result.getDepth(); ++k)
                tmp += ((result.data(i,j,k) >= rgb[k] - threshold[k]) and (result.data(i,j,k) <= rgb[k] + threshold[k]));

            if (tmp == result.getDepth())
                for (int k = 0; k < result.data.depth(); ++k)
                    result.data(i,j,k) = bkg.data(i,j,k);
        }
            
    return result;
}
 
DAISGram DAISGram::equalize() //chiamare grayscale da main
{
    int depth = getDepth(), rows = getRows(), cols = getCols();
    DAISGram result{*this};

    for (int k = 0; k < depth; ++k)
    {
        float vector[256] = {};
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                vector[(int)(*this).data(i,j,k)]++;

        float val{0}, c = 1, cdfmin{};
        for (int i = 0; i < 256; ++i)
            if (vector[i])
            {
                vector[i] += val;
                val = vector[i];
                if(c) cdfmin = vector[i] + (--c);
            }
        
        float den = rows*cols - cdfmin;
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.data(i,j,k) = roundf((((vector[(int)(*this).data(i,j,k)]) - cdfmin) / den) * 255); // v = (int)(*this).data(i,j,k)
    }
    
    return result;
}
 
void DAISGram::generate_random(int h, int w, int d)
{
    data = Tensor(h,w,d,0.0);
    data.init_random(128,50);
    data.rescale(255);
}
