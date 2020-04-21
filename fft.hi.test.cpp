/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/
#include <armadillo>
#include <iostream>
using namespace std;
using namespace arma;
#include "wave2D.h"
#include "my_armadillo.h"

template<typename T1>
float* hilbert1D(T1 *s, int n, float dt);

int main ()
{
   float dt(0.0005),da;
   int i,j, T(512);
   float *s, *sh;
   s=new float[int(T)];
   wavelet01(s, T, dt);
   sh=hilbert1D(s, T, dt);

   ofstream outfs,outfsh;
   outfs.open("signal.bin");
   outfsh.open("sh.bin");
   for(j=0;j<T;j++)
   {
      da=j;
      outfs.write((char*)&s[j],sizeof(float));
      outfs.write((char*)&da,sizeof(float));
      outfsh.write((char*)&sh[j],sizeof(float));
      outfsh.write((char*)&da,sizeof(float));
   }
   outfs.close(); 
   outfsh.close(); 

   fmat as(1,T),ash(1,T),ze(1,T);
   cx_fmat fs(1,T),sfft(1,T);

   ze.fill(0.0);
   for(i=0;i<T;i++)
   {
      as(0,i)=s[i];
      ash(0,i)=sh[i];
   }

   fmat ze2d(T,T);
   cx_fmat s2d(T,T),sfft2d(T,T);
   ze2d.fill(0.0);
   s2d.set_real(ze2d);
   s2d.set_imag(ze2d);
   for(i=0;i<T;i++)
   {
      for(j=0;j<T;j++)
      {
         if((i+j)<T)
         {
            s2d(i,j+i).real(as(0,j));
            //cout<<i+j<<endl;
         }
      }
   }
   for(j=0;j<512;j++)
   {
      for(i=0;i<int(250+0.5*j);i++)
      {
         if(T/2-i>=0)
         s2d(int(250+0.5*j)-i-1,j).real(real(s2d(int(250+0.5*j)-i-1,j))+as(0,T/2-i));
      }
   }
   sfft2d=fft2(s2d,T,T);
   ze2d=real(sfft2d);
   datawrite(ze2d,T,T,"sfft2d.bin");
   ze2d=real(s2d);
   datawrite(ze2d,T,T,"s2d.real.bin");

   for(i=0;i<T;i++)
   {
      for(j=0;j<T;j++)
      {
         if((i+j)<T)
         {
            s2d(i,j+i).imag(ash(0,j));
            //cout<<i+j<<endl;
         }
      }
   }
   for(j=0;j<512;j++)
   {
      for(i=0;i<int(250+0.5*j);i++)
      {
         if(T/2-i>=0)
         s2d(int(250+0.5*j)-i-1,j).imag(imag(s2d(int(250+0.5*j)-i-1,j))+ash(0,T/2-i));
      }
   }
   sfft2d=fft2(s2d,T,T);
   ze2d=real(sfft2d);
   datawrite(ze2d,T,T,"shfft2d.bin");
   ze2d=imag(s2d);
   datawrite(ze2d,T,T,"sh2d.imag.bin");

   for(i=T/2;i<T;i++)
   {
      for(j=0;j<T;j++)
      {
         sfft2d(i,j).real(0.0);
         sfft2d(i,j).imag(0.0);
      }
   }
   cx_fmat sifft2d(T,T);
   sifft2d=ifft2(sfft2d,T,T);
   ze2d=real(sifft2d);
   datawrite(ze2d,T,T,"shifft2d.real.bin");

/*
   fs.set_real(as);
   fs.set_imag(ze);
   sfft.row(0)=fft(fs.row(0),T);
   outfs.open("s.fft.bin");
   for(j=0;j<T;j++)
   {
      da=real(sfft(0,j));
      outfs.write((char*)&da,sizeof(float));
      da=j;
      outfs.write((char*)&da,sizeof(float));
   }
   outfs.close(); 
   
   fs.set_real(as);
   fs.set_imag(ash);
   sfft.row(0)=fft(fs.row(0),T);
   outfs.open("sh.fft.bin");
   for(j=0;j<T;j++)
   {
      da=real(sfft(0,j));
      outfs.write((char*)&da,sizeof(float));
      da=j;
      outfs.write((char*)&da,sizeof(float));
   }
   outfsh.close(); 
*/

   return 0;
}


