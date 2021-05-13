#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>
#include <sstream>

#include "dais_exc.h"
#include "tensor.h"

#define PI 3.141592654
#define FLT_MAX 3.402823466e+38F /* max value */
#define FLT_MIN 1.175494351e-38F /* min positive value */

using namespace std;

Tensor::Tensor()
{
    channels = nullptr;
    r = 0;
    c = 0;
    d = 0;
}

Tensor::Tensor(int r, int c, int d, float v)
{
    if (r <= 0 or c <= 0 or d <= 0) throw dimension_mismatch();

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

Tensor::Tensor(const Tensor& src)
{
    if (src.channels)
    {
        this->r = src.r;
        this->c = src.c;
        this->d = src.d;

        channels = new float*[d];
        channels[0] = new float[r*c*d];

        for (int i = 0; i < d; ++i)
            channels[i] = &channels[0][i*(r*c)];

        for (int i = 0; i < r*c*d; ++i)
            channels[0][i] = src.channels[0][i];
    }
    else
        channels = nullptr;
}

void Tensor::init_random(float mean, float std)
{
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

void Tensor::init(int r, int c, int d, float v)
{
    if (channels)
    {
        delete[] channels[0];
        delete[] channels;
        channels = nullptr;
    }
    (*this) = Tensor(r, c, d, v);
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
    if (not channels) throw tensor_not_initialized(); //caso i,j,k su tensor vuoto
    if (i < 0 or i > r or j < 0 or j > c or k < 0 or k > d)
        throw index_out_of_bound();
    return channels[k][i*c + j];
}

float& Tensor::operator()(int i, int j, int k)
{
    if (not channels) throw tensor_not_initialized();
    if (i < 0 or i > r or j < 0 or j > c or k < 0 or k > d)
        throw index_out_of_bound();
    return channels[k][i*c + j]; //formula array C-style, indice riga * nr. colonne + indice colonna
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
                    throw unknown_operation();
            }

    return result;
}

Tensor Tensor::operator-(const float &rhs) const
{
    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) -= rhs;

    return result;
}

Tensor Tensor::operator+(const float &rhs) const
{
    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) += rhs;

    return result;
}

Tensor Tensor::operator*(const float &rhs) const
{
    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) *= rhs;

    return result;
}

Tensor Tensor::operator/(const float &rhs) const
{
    if (not rhs) throw unknown_operation();

    Tensor result{*this};

    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                result(i,j,k) /= rhs;

    return result;
}

Tensor& Tensor::operator=(const Tensor& other)
{
    if(this == &other) return *this;
    
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

void Tensor::clamp(float low, float high)
{
    for (int k = 0; k < d; ++k)
        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
            {
                if ((*this)(i,j,k) < low) (*this)(i,j,k) = low;
                else if ((*this)(i,j,k) > high) (*this)(i,j,k) = high;
            }
}

void Tensor::rescale(float new_max)
{
    for (int k = 0; k < d; ++k)
    {
        float max = getMax(k);
        float min = getMin(k);

        for (int i = 0; i < r; ++i)
            for (int j = 0; j < c; ++j)
                (*this)(i,j,k) = (((*this)(i,j,k) - min) / (max - min)) * new_max;
    }
}

Tensor Tensor::padding(int pad_h, int pad_w) const
{
    //Chiamata constructor con la nuova dimensione
    Tensor result{r+2*pad_h,c+2*pad_w,d};

    for (int k = 0; k < result.d; ++k) //k for depth, i for rows and j for columns
        for (int i = pad_h; i < result.r - pad_h; ++i)
            for (int j = pad_w; j < result.c - pad_w; ++j)
                result(i,j,k) = (*this)(i-pad_h,j-pad_w,k);

    return result;
}

Tensor Tensor::subset(int row_start, 
                      int row_end, 
                      int col_start, 
                      int col_end, 
                      int depth_start, 
                      int depth_end) const
{
    if (not channels) throw tensor_not_initialized();
    if (row_start > r or row_end > r or row_start < 0 or row_end < 0 or
        col_start > c or col_end > c or col_start < 0 or col_end < 0 or
        depth_start > d or depth_end > d or depth_start < 0 or depth_end < 0)
            throw dimension_mismatch();
    
    Tensor result{row_end - row_start, col_end - col_start, depth_end - depth_start};

    for (int k = depth_start; k < depth_end; ++k)
        for (int i = row_start; i < row_end; ++i)
            for (int j = col_start; j < col_end; ++j)
                result(i - row_start, j - col_start, k - depth_start) = (*this)(i,j,k);

    return result;
}

Tensor Tensor::concat(const Tensor &rhs, int axis) const 
{
    Tensor result{};
    switch (axis)
    {
        case 0: (c == rhs.c and d == rhs.d) ? result.init(r + rhs.r, c, d) : throw concat_wrong_dimension(); break;
        case 1: (r == rhs.r and d == rhs.d) ? result.init(r, c + rhs.c, d) : throw concat_wrong_dimension(); break;
        case 2: (c == rhs.c and r == rhs.r) ? result.init(r, c, d + rhs.d) : throw concat_wrong_dimension(); break;
        default: throw concat_wrong_dimension();
    }

    for (int k = 0; k < result.d; ++k) //k for depth, i for rows and j for columns
        for (int i = 0; i < result.r; ++i)
            for (int j = 0; j < result.c; ++j)
            {
                switch (axis)
                {
                    case 1: (j >= c) ? result(i,j,k) = rhs(i,j-c,k) : result(i,j,k) = (*this)(i,j,k); break;
                    case 2: (k >= d) ? result(i,j,k) = rhs(i,j,k-d) : result(i,j,k) = (*this)(i,j,k); break;
                    default: (i >= r) ? result(i,j,k) = rhs(i-r,j,k) : result(i,j,k) = (*this)(i,j,k); break; //default case is axis == 0
                }
            }

    return result;
}

Tensor Tensor::convolve(const Tensor &f) const
{
    if (f.r % 2 == 0 or f.c % 2 == 0) throw filter_odd_dimensions();

    int pad_h = (f.r - 1) / 2;
    int pad_w = (f.c - 1) / 2;
    Tensor result{r,c,d}, tmp_tensor{};
    int col_start{0};
    float tmp_float{0};

    tmp_tensor = padding(pad_h,pad_w);

    for (int i = 0; i < r; ++i)
    {
        for (int j = 0; j < c; ++j)
        {
            for (int k = 0; k < d; k++)
            {
                for (int i2 = i; i2 < i + f.r; ++i2)
                    for (int j2 = col_start; j2 < col_start + f.c; ++j2)
                        tmp_float += tmp_tensor(i2,j2,k) * f(i2 - i, j2 - col_start, 0);

                result(i,j,k) = tmp_float;
                tmp_float = 0;
            }
            col_start++;
        }
        col_start = 0;
    }

    return result;
}

int Tensor::rows()const {return r;}

int Tensor::cols()const {return c;}

int Tensor::depth()const {return d;}

float Tensor::getMin(int k) const
{
    if (not channels) throw tensor_not_initialized();
    if (k > d or k < 0) throw index_out_of_bound();

    float min{(*this)(0,0,k)};

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if ((*this)(i,j,k) < min) min = (*this)(i,j,k);
        
    return min;
}

float Tensor::getMax(int k) const
{
    if (not channels) throw tensor_not_initialized();
    if (k > d or k < 0) throw index_out_of_bound();

    float max{(*this)(0,0,k)};

    for (int i = 0; i < r; ++i)
        for (int j = 0; j < c; ++j)
            if ((*this)(i,j,k) > max) max = (*this)(i,j,k);
        
    return max;
}

void Tensor::showSize()const
{
    cout << r << " x " << c << " x " << d << endl;
}

ostream& operator<<(ostream& stream, const Tensor & obj)
{
    if(not obj.channels) throw tensor_not_initialized();

    for(int k = 0; k < obj.d; k++)
    {
        for(int i = 0; i < obj.r; i++)
        {
            for(int j = 0; j < obj.c; j++)
                stream << "[" << obj(i,j,k) << "] ";
            stream << endl;
        }
        stream << endl;
    }
    return stream;
}

void Tensor::read_file(string filename) // dubbio: nella write_file inserisce le quadre, ma nel read non le considera.  non rischia di dare problemi?
{
    ifstream myfile(filename); //open file

    if (not myfile.is_open()) throw unable_to_read_file();
    
    myfile >> r >> c >> d; //prime tre righe

    init(r, c, d);

    for(int k = 0; k < d; k++)
        for(int i = 0; i < r; i++)
            for(int j = 0; j < c; j++)
                myfile >> (*this)(i,j,k); 
    
    myfile.close();
}

void Tensor::write_file(string filename)
{
    ofstream myfile(filename); //open file

    myfile << r << endl << c << endl << d << endl; //prime tre righe

    for(int k = 0; k < d; k++)
        for(int i = 0; i < r; i++)
            for(int j = 0; j < c; j++)
                myfile << (*this)(i,j,k) << endl;
    
    myfile.close();
}