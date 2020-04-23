/*********(version 1.0)***********/
/*注释�?
    C++程序模板�?

*/
/********************************/

#include <iostream>
using namespace std;

#include <xj.c++.h>

float Blackman(int n, int N)
{
    float xs, pi(3.1415926);
    xs=0.42-0.5*cos(2*pi*n/N/2)+0.08*cos(4*pi*n/N/2);
    return xs;
}

cx_fmat windowsZ(cx_fmat &sfft, int Z, int X);
void windowsSource(fmat &cp, int Z, int X, int zf)
{
    int i,j,n,N(10);
    float xs;
    for(i=0;i<Z;i++)
    {
        for(j=0;j<X;j++)
        {
            if(i<zf && zf-i<=N)
            {
                xs=Blackman((zf-i),N);
                cp(i,j)=cp(i,j)*xs;
            }
            else if(i>zf && i-zf<=N)
            {
                xs=Blackman((i-zf),N);
                cp(i,j)=cp(i,j)*xs;
            }
            else if(i==zf)
            {
                cp(i,j)=0.0;
            }
        }
    }
}

int main ()
{
    int i,j,k,T(2000),Z(200),Zfft(256),X(401),Xfft(512),sx(200),ds(10),s_id,s_num(1),zf(35);
    float *s, *sh;
    fmat suf(T,X),suf_v(T,X),a(Z,X),b(Z,X),c(Z,X),d(Z,X);
    char file1[99],file2[99],file3[99],file4[99],file5[99];
    fmat sufh(T,X),sufh_v(T,X),ns(T,X),nsh(T,X);

    ofstream outfs, outfsh, outfup, outfdown;
    wave2D S(Z,X), Sh(Z,X);
    fmat cp(Z,X),cpfft(Zfft,Xfft);
    cx_fmat comwave(Z,X), wfft(Zfft,Xfft), wfft2(Zfft,Xfft), wifft(Zfft,Xfft);

    s=new float[T];
    wavelet01(s, T, S.dt);
    sh=hilbert1D(s, T, S.dt);
    S.cleardata();
    Sh.cleardata();
    dataread(S.p2,Z,X,"./data/model.dat");  
    dataread(Sh.p2,Z,X,"./data/model.dat");  

    outfs.open("./data/zy.movie.s.bin");
    outfsh.open("./data/zy.movie.sh.bin");
    outfup.open("./data/zy.movie.up.bin");
    outfdown.open("./data/zy.movie.down.bin");
    for(k=0;k<T;k++)
        {
        S.s2[zf][sx]+=s[k];
        S.timeslicecal();
        S.timeslicecopy();
        Sh.s2[zf][sx]+=sh[k];
        Sh.timeslicecal();
        Sh.timeslicecopy();
        
        datawrite(S.s3, Z, X, outfs);
        datawrite(Sh.s3, Z, X, outfsh);

        cp=matcopy(S.s3,Z,X);
        //windowsSource(cp, Z, X, zf);
        comwave.set_real(cp);
        cp=matcopy(Sh.s3,Z,X);
        //windowsSource(cp, Z, X, zf);
        comwave.set_imag(cp);
        wfft=fft2(comwave,Zfft,Xfft);
        //wfft=windowsZ(wfft, Z, X);

        cpfft.fill(0.0);
        wfft2.set_real(cpfft);
        wfft2.set_imag(cpfft);
        for(i=0;i<Zfft/2;i++)
        {
            wfft2.row(i)=wfft.row(i);
        }
        wifft=ifft2(wfft2,Zfft,Xfft);
        cpfft=real(wifft);
        datawrite(cpfft,Zfft,Xfft,outfup);

        cpfft.fill(0.0);
        wfft2.set_real(cpfft);
        wfft2.set_imag(cpfft);
        for(i=Zfft/2;i<Zfft;i++)
        {
            wfft2.row(i)=wfft.row(i);
        }
        wifft=ifft2(wfft2,Zfft,Xfft);
        cpfft=real(wifft);
        datawrite(cpfft,Zfft,Xfft,outfdown);

        b=matcopy(S.s3,Z,X);
        suf.row(k)=b.row(zf);
        b=matcopy(Sh.s3,Z,X);
        sufh.row(k)=b.row(zf);
        if(k%100==0)
            cout<<k<<endl;
        }
    outfs.close();
    outfsh.close();
    outfup.close();
    outfdown.close();

    file2[0]='\0';
    strcat(file2,"./data/sufs.bin");
    //strcat(file2,numtostr(s_id,5));
    datawrite(suf,T,X,file2);

    file2[0]='\0';
    strcat(file2,"./data/sufsh.bin");
    //strcat(file2,numtostr(s_id,5));
    datawrite(sufh,T,X,file2);

    S.cleardata();
    S.setvelocity(3000);
    Sh.cleardata();
    //wave2D C(Z,X);
    Sh.setvelocity(3000);
    for(k=0;k<T;k++)
    {
        S.s2[zf][sx]+=s[k];
        S.timeslicecal();
        S.timeslicecopy();
        Sh.s2[zf][sx]+=sh[k];
        Sh.timeslicecal();
        Sh.timeslicecopy();

        b=matcopy(S.s3,Z,X);
        suf_v.row(k)=b.row(zf);
        b=matcopy(Sh.s3,Z,X);
        sufh_v.row(k)=b.row(zf);
        if(k%100==0)
            cout<<k<<endl;
    }

    file3[0]='\0';
    strcat(file3,"./data/sufsv.bin");
    //strcat(file3,numtostr(s_id,5));
    datawrite(suf_v,T,X,file3);
    
    file3[0]='\0';
    strcat(file3,"./data/sufshv.bin");
    //strcat(file3,numtostr(s_id,5));
    datawrite(sufh_v,T,X,file3);

    S.cleardata();
    Sh.cleardata();
    dataread(S.p2,Z,X,"./data/model.dat");  
    dataread(Sh.p2,Z,X,"./data/model.dat");  
    dataread(suf_v,T,X,"./data/sufsv.bin");
    dataread(ns,T,X,"./data/sufs.bin"); 
    ns=ns-suf_v;
    dataread(suf_v,T,X,"./data/sufshv.bin");
    dataread(nsh,T,X,"./data/sufsh.bin"); 
    nsh=nsh-suf_v;

    outfs.open("./data/ns.movie.s.bin");
    outfsh.open("./data/ns.movie.sh.bin");
    outfup.open("./data/ns.movie.up.bin");
    outfdown.open("./data/ns.movie.down.bin");
    for(k=0;k<T;k++)
        {
        d=matcopy(S.s2,Z,X);
        d.row(zf)=d.row(zf)+ns.row(T-1-k);
        matcopy(S.s2,d,Z,X);
        S.timeslicecal();
        S.timeslicecopy();
        d=matcopy(Sh.s2,Z,X);
        d.row(zf)=d.row(zf)+nsh.row(T-1-k);
        matcopy(Sh.s2,d,Z,X);
        Sh.timeslicecal();
        Sh.timeslicecopy();
        datawrite(S.s3, Z, X, outfs);
        datawrite(Sh.s3, Z, X, outfsh);

        cp=matcopy(S.s3,Z,X);
        //windowsSource(cp, Z, X, zf);
        comwave.set_real(cp);
        cp=matcopy(Sh.s3,Z,X);
        //windowsSource(cp, Z, X, zf);
        comwave.set_imag(cp);
        wfft=fft2(comwave,Zfft,Xfft);
        //wfft=windowsZ(wfft, Z, X);

        cpfft.fill(0.0);
        wfft2.set_real(cpfft);
        wfft2.set_imag(cpfft);
        for(i=0;i<Zfft/2;i++)
        {
            wfft2.row(i)=wfft.row(i);
        }
        wifft=ifft2(wfft2,Zfft,Xfft);
        cpfft=real(wifft);
        datawrite(cpfft,Zfft,Xfft,outfup);

        cpfft.fill(0.0);
        wfft2.set_real(cpfft);
        wfft2.set_imag(cpfft);
        for(i=Zfft/2;i<Zfft;i++)
        {
            wfft2.row(i)=wfft.row(i);
        }
        wifft=ifft2(wfft2,Zfft,Xfft);
        cpfft=real(wifft);
        datawrite(cpfft,Zfft,Xfft,outfdown);
        
        if(k%100==0)
            cout<<k<<endl;
        }
    outfs.close();
    outfsh.close();
    outfup.close();
    outfdown.close();

    return 0;
}

cx_fmat windowsZ(cx_fmat &sfft, int Z, int X)
{
    int i, j, N;
    float k(0.5), xs;
    
    for(i=0;i<Z;i++)
    {
        for(j=0;j<X;j++)
        {
            if(i<Z/2)
            {   
                N=int(i*k)+1;
                if(j<X/2 && j<=N)
                    {
                        xs=Blackman(j,N);
                        sfft(i,j).real(real(sfft(i,j))*xs);
                        sfft(i,j).imag(imag(sfft(i,j))*xs);
                    }
                if(j>=X/2 && X-1-j<=N)
                    {
                        xs=Blackman(X-1-j,N);
                        sfft(i,j).real(real(sfft(i,j))*xs);
                        sfft(i,j).imag(imag(sfft(i,j))*xs);
                    }
            }
            if(i>=Z/2)
            {   
                N=int((Z-1-i)*k)+1;
                if(j<X/2 && j<=N)
                    {
                        xs=Blackman(j,N);
                        sfft(i,j).real(real(sfft(i,j))*xs);
                        sfft(i,j).imag(imag(sfft(i,j))*xs);
                    }
                if(j>=X/2 && X-1-j<=N)
                    {
                        xs=Blackman(X-1-j,N);
                        sfft(i,j).real(real(sfft(i,j))*xs);
                        sfft(i,j).imag(imag(sfft(i,j))*xs);
                    }
            }

        }
    }
    return sfft;
}


