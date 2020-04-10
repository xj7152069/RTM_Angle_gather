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
    float f,**m;
    fmat suf(T,X),suf_v(T,X),a(Z,X),b(Z,X),c(Z,X),ns(T,X),d(Z,X);
    char file1[99],file2[99],file3[99],file4[99],file5[99];

    ofstream outf1;
    wave2D A(Z,X);
    
    for(s_id=0;s_id<s_num;s_id++)
    {
        cout<<sx<<" : "<<endl;

        A.cleardata();
        dataread(A.p2,Z,X,"model.dat");  
        //dataread(m,Z,X,"model.dat");
        //matmul(m,0.000000001,Z,X); 
        //a=matcopy(m,Z,X);
        file1[0]='\0';
        strcat(file1,"movie");
        strcat(file1,numtostr(s_id,5));
        outf1.open(file1);
        for(k=0;k<T;k++)
            {
            f=wavelet01(k,A.dt);
            A.s2[zf][sx]=A.s2[zf][sx]+f;
            A.timeslicecal();
            A.timeslicecopy();
            
            datawrite(A.s3, Z, X, outf1);
            b=matcopy(A.s3,Z,X);
            suf.row(k)=b.row(zf);
            }
        outf1.close();

        file2[0]='\0';
        strcat(file2,"suf");
        strcat(file2,numtostr(s_id,5));
        datawrite(suf,T,X,file2);

        A.cleardata();
        //wave2D C(Z,X);
        A.setvelocity(3000);
        for(k=0;k<T;k++)
            {
            f=wavelet01(k,A.dt);
            A.s2[zf][sx]=A.s2[zf][sx]+f;
            A.timeslicecal();
            A.timeslicecopy();
            
            b=matcopy(A.s3,Z,X);
            suf.row(k)=b.row(zf);
            }

        file3[0]='\0';
        strcat(file3,"sufv");
        strcat(file3,numtostr(s_id,5));
        datawrite(suf,T,X,file3);

        A.cleardata();
        //wave2D B(Z,X);
        dataread(A.p2,Z,X,"model.dat");  
        dataread(suf_v,T,X,file3);
        dataread(ns,T,X,file2); 
        ns=ns-suf_v;

        file4[0]='\0';
        strcat(file4,"nspy");
        strcat(file4,numtostr(s_id,5));
        datawrite(ns,T,X,file4);

        file5[0]='\0';
        strcat(file5,"nsmovie");
        strcat(file5,numtostr(s_id,5));
        outf1.open(file5);
        for(k=0;k<T;k++)
            {
            d=matcopy(A.s2,Z,X);
            d.row(zf)=d.row(zf)+ns.row(T-1-k);
            matcopy(A.s2,d,Z,X);
            A.timeslicecal();
            A.timeslicecopy();
            
            datawrite(A.s3, Z, X, outf1);
            }
        outf1.close();
        //B.~wave2D();
        sx=sx+ds;
    }

    return 0;
}


