/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

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
    Poynting_2D A(3,nz,nx);

    ofstream out1,out2,out3;
    ifstream inf1;

    inf1.open("nsmovie0.000");
    out1.open("vx.ns.bin");
    out2.open("vz.ns.bin");
    out3.open("wave.ns.bin");
    
    for(k=1620;k<1625;k++)
    {
        inf1.seekg(nz*nx*k*sizeof(float), ios::beg);
        dataread(p,nz,nx,inf1);
        
        A.addtimeslicecal(p);
        
        A.velocityCalculate();
        if(k==1623)
        {
            matsmooth(vx,A.vx,nz,nx);
            matsmooth(vz,A.vz,nz,nx);
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
