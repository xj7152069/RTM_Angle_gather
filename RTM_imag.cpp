/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
using namespace std;

#include "wave2D.h"
#include "my_armadillo.h"

int main ()
{
   int k,T(2000),Z(200),X(401),sx(200),ds(10),s_id,s_num(1),zf(35);
   fmat imag(Z,X),imagla(Z,X),t1(Z,X),t2(Z,X);
   fmat la(Z,X);
   char file1[99],file2[99],file3[99],file4[99],file5[99];

   ifstream inf,inf_ns;
   imag.fill(0.0);
   imagla.fill(0.0);

   int i,j;

   for(s_id=0;s_id<s_num;s_id++)
   {
      cout<<s_id<<" : "<<endl;
      file1[0]='\0';
      strcat(file1,"movie");
      strcat(file1,numtostr(s_id,5));
      inf.open(file1);

      file5[0]='\0';
      strcat(file5,"nsmovie");
      strcat(file5,numtostr(s_id,5));
      inf_ns.open(file5);

      for(k=0;k<T;k++)
      {
         inf_ns.seekg((T-k-1)*Z*X*sizeof(float),ios::beg);
         dataread(t1,Z,X,inf);
         dataread(t2,Z,X,inf_ns);
         t2=matmul(t2,t1,Z,X);
         imag=imag+t2;
      }
      inf.close();
      inf_ns.close();
   }

   for(i=1;i<Z-1;i++)
      {
         for(j=1;j<X-1;j++)
         {
            imagla(i,j)=(imag(i-1,j-1)-imag(i,j)+imag(i-1,j+1)\
            -imag(i,j)+imag(i+1,j-1)-imag(i,j)+imag(i+1,j+1)-imag(i,j));
         }
      }
   datawrite(imag,Z,X,"imag.bin");
   datawrite(imagla,Z,X,"imag2.bin");

   return 0;
}




