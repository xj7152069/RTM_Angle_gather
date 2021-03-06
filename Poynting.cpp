/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
using namespace std;
 
#include <xj.c++.h>

int main ()
{
    int nz(200),nx(401),nt(2000),na(200),xfft(512),zfft(256);
    int i,j,k,n,N(25);
    float **vx,**vz,**swave,**rwave,**imag;
    swave=newfmat(zfft,xfft);
    rwave=newfmat(zfft,xfft);
    imag=newfmat(nz,nx);
    vx=newfmat(nz,nx);
    vz=newfmat(nz,nx);
    //matcopy(p,0.0,nz,nx);
    matcopy(imag,0.0,nz,nx);
    matcopy(vx,0.0,nz,nx);
    matcopy(vz,0.0,nz,nx);
    Poynting_2D S(3,nz,nx);
    Poynting_2D R(3,nz,nx);
    angle_gather2D test(na, nz, nx);

    //ofstream out1,out2,out3;
    ifstream infs,infr;
    ofstream outfcig,outfimag;
    ofstream outstheta,outrtheta;

    infs.open("./data/zy.movie.down.bin");
    infr.open("./data/ns.movie.up.bin");
    outstheta.open("./data/thetaS.down.bin");
    outrtheta.open("./data/thetaR.up.bin");
    //out3.open("./data/wave.ns.pyt");
    
    for(k=0;k<2000;k++)
    {
        infr.seekg(zfft*xfft*(2000-k-1)*sizeof(float), ios::beg);
        dataread(rwave,zfft,xfft,infr);
        dataread(swave,zfft,xfft,infs);
        
        S.addtimeslicecal(swave);
        R.addtimeslicecal(rwave);
        
        S.velocityCalculate();
        R.velocityCalculate();

        matsmooth(vz,S.vz,nz,nx);
        matsmooth(vx,S.vx,nz,nx);
        test.theatcal_S(vz, vx);
        matsmooth(vz,R.vz,nz,nx);
        matsmooth(vx,R.vx,nz,nx);
        test.theatcal_R(vz, vx);
        
        datawrite(test.pr,nz,nx,outrtheta);
        datawrite(test.ps,nz,nx,outstheta);

        test.addFA(swave, rwave);

        /*
        if(k==1623)
        {
            //matsmooth(vx,A.vx,nz,nx);
            //matsmooth(vz,A.vz,nz,nx);
            //velocitynormal(vx,vz,nz,nx);
            datawrite(vx,nz,nx,out1);
            datawrite(vz,nz,nx,out2);
            datawrite(p,nz,nx,out3);
        }*/
        
        if(k%100==0)
            cout<<k<<endl;
    }
    //out2.close();
    //out1.close();
    outfcig.open("./data/fa.angle.cig.gather");
    for(i=0;i<nx;i++)
    {
        datawrite(test.FA[i],nz,na,outfcig);
        //cout<<"output the CIG"<<endl;
        for(j=0;j<nz;j++)
        {
            for(k=0;k<na;k++)
            {
                imag[j][i]+=test.FA[i][j][k];
            }
        }
    }
    datawrite(imag,nz,nx,"./data/fa.angle.imag");

   return 0;
}
