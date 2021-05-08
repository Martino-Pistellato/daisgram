#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

void Tensor::init_random(float mean, float std){
    if(channels){

        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i=0;i<r;i++){
            for(int j=0;j<c;j++){
                for(int k=0;k<d;k++){
                    this->operator()(i,j,k)= distribution(generator);
                }
            }
        }    

    }else{
        throw(tensor_not_initialized());
    }
}

Tensor::Tensor(int r, int c, int d, float v)
{
    this->r = r;
    this->c = c;
    this->d = d;
    channels = new float**[d];
    channels[0] = new float*[r*d];
    channels[0][0] = new float[r*c*d]; //matrix is the big flattened array
    
    for (int i = 0; i < d; ++i)
    {
        channels[i] = &channels[0][i*r];

        for (int j = 0; j < r; ++j)
        {
            channels[0][j+(i*r)] = &channels[0][0][(j+(i*r))*c];
            for (int k = 0; k < c; k++)
                channels[0][0][(j+(i*r))*c + k] = v;
        }
    }
}

Tensor::Tensor(const Tensor& that)
{
    delete[] channels[0][0];
    delete[] channels[0];
    delete[] channels;

    if (that.channels)
    {
        this->r = that.r;
        this->c = that.c;
        this->d = that.d;

        channels = new float**[d];
        channels[0] = new float*[r*d];
        channels[0][0] = new float[r*c*d];

        for (int i = 0; i < d; ++i)
            for (int j = 0; j < r; ++j)
                for (int k = 0; k < c; ++k)
                    channels[i][j][k] = that.channels[i][j][k];
    }
    else
        channels = nullptr;
}

Tensor::~Tensor()
{
    delete[] channels[0][0];
    delete[] channels[0];
    delete[] channels;
}

float Tensor::operator()(int i, int j, int k) const
{
    if (i < 0 or i > r or j < 0 or j > c or k < 0 or k > d)
        throw index_out_of_bound();
    return channels[k][i][j];
}

float& Tensor::operator()(int i, int j, int k)
{
    if (i < 0 or i > r or j < 0 or j > c or k < 0 or k > d)
        throw index_out_of_bound();
    return channels[k][i][j];
}

int Tensor::rows() {return r;}

int Tensor::depth() {return d;}

float Tensor::getMax(int k)
{
    if (not channels) throw tensor_not_initialized();
    if (k > d or k < 0) throw index_out_of_bound();

    float max{channels[k][0][0]};

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if (channels[k][i][j] > max) max = channels[k][i][j];
    
    return max;
}

float Tensor::getMin(int k)
{
    if (not channels) throw tensor_not_initialized();
    if (k > d or k < 0) throw index_out_of_bound();

    float min{channels[k][0][0]};

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if (channels[k][i][j] < min) min = channels[k][i][j];
    
    return min;
}

Tensor Tensor::operator-(const Tensor &rhs)
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < r; ++j)
            for (int k = 0; k < c; ++k)
                result(j,k,i) -= rhs(j,k,i);
    
    return result;
}

Tensor Tensor::operator+(const Tensor &rhs)
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < r; ++j)
            for (int k = 0; k < c; ++k)
                result(j,k,i) += rhs(j,k,i);
    
    return result;
}

Tensor Tensor::operator*(const Tensor &rhs)
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < r; ++j)
            for (int k = 0; k < c; ++k)
                result(j,k,i) *= rhs(j,k,i);
    
    return result;
}

Tensor Tensor::operator/(const Tensor &rhs)
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int i = 0; i < d; ++i)
        for (int j = 0; j < r; ++j)
            for (int k = 0; k < c; ++k)
            {
                if (rhs(j,k,i))
                    result(j,k,i) /= rhs(j,k,i);
                else
                    throw unknown_operation();
            }
    
    return result;
}