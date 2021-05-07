#include <iostream>
#include <string>
#include <random>
#include <math.h>
#include <fstream>

using namespace std;

class Tensor
{
    private:

        float *** channels;

        int c;
        int r; 
        int d;

    public:
    
    Tensor(int r, int c, int d, float v = 0.0)
    {
        this->r = r;
        this->c = c;
        this->d = d;

        float*** channels = new float**[d];
        channels[0] = new float*[r*d];
        channels[0][0] = new float[r*c*d]; //it is the big flattened array
        
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

    ~Tensor()
    {
        delete[] channels[0][0];
        delete[] channels[0];
        delete[] channels;
    }
        
    void stampa()
    {
        for (int i = 0; i < d; i++)
        {
            for (int j = 0; j < r; j++)
            {
                for(int k = 0; k < c; k++)
                    cout << channels[i][j][k] << " ";
                cout << endl;
            }
            cout << endl;
        } 
    }
        
        
};

int main()
{
    Tensor prova(3,3,3);
    
    prova.stampa();
    return 0;
}