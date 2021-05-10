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

void Tensor::init_random(float mean, float std)
{
    if (channels)
    {
        std::default_random_engine generator;
        std::normal_distribution<float> distribution(mean,std);

        for(int i = 0; i < r; ++i)
            for(int j = 0; j < c; ++j)
                for(int k = 0; k < d; ++k)
                    this->operator()(i,j,k) = distribution(generator);
    }
    else
        throw(tensor_not_initialized());
}

Tensor::Tensor()
{
    channels = nullptr;
    r = 0;
    c = 0;
    d = 0;
}

Tensor::Tensor(int r, int c, int d, float v)
{
    if (r <= 0 or c <= 0 or d <= 0)
            throw dimension_mismatch();

        this->r = r;
        this->c = c;
        this->d = d;

        channels = new float*[d];
        channels[0] = new float[r*c*d];

        for (int i = 0; i < d; ++i)
            channels[i] = &channels[0][i*(r*c)];

        for (int i = 0; i < r*c*d; ++i)
            channels[0][i] = v;
}

Tensor::Tensor(const Tensor& that)
{
    if (that.channels)
    {
        this->r = that.r;
        this->c = that.c;
        this->d = that.d;

        channels = new float*[d];
        channels[0] = new float[r*c*d];

        for (int i = 0; i < d; ++i)
            channels[i] = &channels[0][i*(r*c)];

        for (int i = 0; i < r*c*d; ++i)
            channels[0][i] = that.channels[0][i];
    }
    else
        channels = nullptr;
}

Tensor::~Tensor()
{
    if (channels)
    {
        delete[] channels[0];
        delete[] channels;
        channels = nullptr;
    }
}

float Tensor::operator()(int i, int j, int k) const
{
    if (i < 0 or i > r or j < 0 or j > c or k < 0 or k > d)
        throw index_out_of_bound();
    return channels[k][i*c + j];
}

float& Tensor::operator()(int i, int j, int k)
{
    if (i < 0 or i > r or j < 0 or j > c or k < 0 or k > d)
        throw index_out_of_bound();
    return channels[k][i*c + j];
}

float Tensor::getMax(int k) const
{
    if (not channels) throw tensor_not_initialized();
    if (k > d or k < 0) throw index_out_of_bound();

    float max{channels[k][0]};

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if (channels[k][i*r + j] > max) max = channels[k][i*r + j];
        
    return max;
}

float Tensor::getMin(int k) const
{
    if (not channels) throw tensor_not_initialized();
    if (k > d or k < 0) throw index_out_of_bound();

    float min{channels[k][0]};

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if (channels[k][i*r + j] < min) min = channels[k][i*r + j];
        
    return min;
}

Tensor Tensor::operator-(const Tensor &rhs) const
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) -= rhs(i,j,k);
        
    return result;
}

Tensor Tensor::operator+(const Tensor &rhs) const
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) += rhs(i,j,k);
        
    return result;
}

Tensor Tensor::operator*(const Tensor &rhs) const
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) *= rhs(i,j,k);
        
    return result;
}

Tensor Tensor::operator/(const Tensor &rhs) const
{
    if (r != rhs.r or c != rhs.c or d != rhs.d) throw dimension_mismatch();

    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
            {
                if (rhs(i,j,k))
                    result(i,j,k) /= rhs(i,j,k);
                else
                    throw unknown_operation{};
            }

    return result;
}

Tensor& Tensor::operator=(const Tensor& other)
{
    if (channels)
    {
        delete[] channels[0];
        delete[] channels;
    }
    if (other.channels)
    {
        this->r = other.r;
        this->c = other.c;
        this->d = other.d;

        channels = new float*[d];
        channels[0] = new float[r*c*d];

        for (int i = 0; i < d; ++i)
            channels[i] = &channels[0][i*(r*c)];

        for (int i = 0; i < r*c*d; ++i)
            channels[0][i] = other.channels[0][i];
    }
    else
        channels = nullptr;

    return *this;
}

void Tensor::init(int r, int c, int d, float v)
{
    Tensor(r, c, d, v);
}

Tensor Tensor::padding(int pad_h, int pad_w)const
{
    //Chiamata constructor con la nuova dimensione
    Tensor result{r+2*pad_h,c+2*pad_w,d};

    for (int k = 0; k < d; ++k) //k for depth, i for rows and j for columns
        for (int i = pad_h; i < result.r - pad_h; ++i)
            for (int j = pad_w; j < result.c - pad_w; ++j)
                result(i,j,k) = (*this)(i-pad_h,j-pad_w,k);

    return result;
}