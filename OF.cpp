#include <iostream>
using namespace std;

#include "wave2D.h"
#include "RTMangle2D.h"

void matsmooth(float **mat1, float **mat2, int x1, int x2)
{
    int i,j;
    for(i=1;i<x1-1;i++)
    {
        for(j=1;j<x2-1;j++)
        {
            mat1[i][j]=(mat2[i-1][j-1]*1.0/12+mat2[i-1][j]*1.0/6+mat2[i-1][j+1]*1.0/12\
                    +mat2[i][j-1]*1.0/6+mat2[i][j+1]*1.0/6+mat2[i+1][j-1]*1.0/12\
                    +mat2[i+1][j]*1.0/6+mat2[i+1][j+1]*1.0/12+mat2[i][j]*1.0/3);
        }
    }
}

void velocitynormal(float **vx, float **vz, int x1, int x2)
{
    int i,j;
    float vl;
    for(i=1;i<x1-1;i++)
    {
        for(j=1;j<x2-1;j++)
        {
            vl=sqrt(vx[i][j]*vx[i][j]+vz[i][j]*vz[i][j]);
            if(vl>0.000000000001)
            {
                vx[i][j]=vx[i][j]/vl;
                vz[i][j]=vz[i][j]/vl;
            }
        }
    }
}

int main ()
{
    int nz(200),nx(401),nt(2000);
    int i,j,k,n,N(25);
    float **vx,**vz,**p;
    p=newfmat(nz,nx);
    vx=newfmat(nz,nx);
    vz=newfmat(nz,nx);
    matcopy(p,0.0,nz,nx);
    matcopy(vx,0.0,nz,nx);
    matcopy(vz,0.0,nz,nx);
    OF_2D A(3,nz,nx);

    ofstream out1,out2,out3;
    ifstream inf1;

    inf1.open("movie0.000");
    out1.open("vx.zy.bin");
    out2.open("vz.zy.bin");
    out3.open("wave.zy.bin");
    
    for(k=1620;k<1625;k++)
    {
        inf1.seekg(nz*nx*k*sizeof(float), ios::beg);
        dataread(p,nz,nx,inf1);
        
        A.addtimeslicecal(p);
        
        A.velocityInitialization();
        //cout<<"ok"<<endl;
        for(n=0;n<N;n++)
        {
            A.velocityIteration();
        }
        if(k==1623)
        {
            matsmooth(vx,A.vx2,nz,nx);
            matsmooth(vz,A.vz2,nz,nx);
            //velocitynormal(vx,vz,nz,nx);
            datawrite(vx,nz,nx,out1);
            datawrite(vz,nz,nx,out2);
            datawrite(p,nz,nx,out3);
        }
        
        if(k%100==0)
            cout<<k<<endl;
    }
    out2.close();
    out1.close();

    return 0;
}



