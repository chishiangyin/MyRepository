//可能需添加错误输入（如负能量的保护）
#include "BetaDecayBranch.hh"
#include <iostream>
#include <cmath>
#include <string>
#include "TMath.h"
#include "TF1.h"
BetaDecayBranch::BetaDecayBranch()
{
    alpha=1.0/137.03604;                  
    me=0.511099895;
    mN=939.5654205;
    mu_v = 4.7;
    ga = 1.27590;
}
BetaDecayBranch::~BetaDecayBranch()
{

}
double BetaDecayBranch::FermiFunc(int Z,int A,double Ee)
{
    double R = NuclearRadius(A);
    double pe=sqrt(Ee*Ee-me*me);
    double gamma=sqrt(1-(alpha*double(Z))*(alpha*double(Z)));
    double y=alpha*double(Z)*Ee/pe;
    double F;
    //std::cout<<Ee<<" "<<R<<" "<<pe<<" "<<gamma<<" "<<y<<" "<<std::endl;
    F=2*(gamma+1)*pow(2*pe*R,2*(gamma-1))*exp(M_PI*y)*BIG_ComplexGammaFuc(gamma,y)/((tgamma(2*gamma+1))*(tgamma(2*gamma+1)));
    return F;
}
//Double beta Decay and Majorana Neutrino
double BetaDecayBranch::FermiFuncN(int Z,int A,double Ee,int n)
{
    double R = NuclearRadius(A);
    double pe=sqrt(Ee*Ee-me*me);
    double gamman=sqrt(double(n+1)*double(n+1)-(alpha*double(Z))*(alpha*double(Z)));
    double y=alpha*double(Z)*Ee/pe;
    double F;
    F=pow(2*pe*R,2*(gamman-double(n+1)))*exp(M_PI*y)*BIG_ComplexGammaFuc(gamman,y)*(tgamma(2*double(n+1))/tgamma(double(n+1))/tgamma(2*y+1))*(tgamma(2*double(n+1))/tgamma(double(n+1))/tgamma(2*y+1));
    return F;
}
//New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
double BetaDecayBranch::P32FermiLikeFunc(int Z,int A,double Ee)
{
    return FermiFuncN(Z,A,Ee,1)/FermiFuncN(Z,A,Ee,0);
}
//New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
double BetaDecayBranch::P12FermiLikeFunc(int Z,int A,double Ee)
{
    double R = NuclearRadius(A);
    double pe = sqrt(Ee*Ee-me*me);
    return ((alpha*Z/2.0+Ee*R/3.0)*(alpha*Z/2.0+Ee*R/3.0)+(me*R/3.0)*(me*R/3.0)-2*me*me*R*(alpha*Z/2.0+Ee*R/3.0)/(3.0*Ee))/(std::sph_bessel(1,pe*R)*std::sph_bessel(1,pe*R));
}
//New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
double BetaDecayBranch::SP12FermiLikeFunc(int Z,int A,double Ee)
{
    double R = NuclearRadius(A);
    double pe = sqrt(Ee*Ee-me*me);
    return ((alpha*Z/2.0+Ee*R/3.0)-(me*me*R/3.0/Ee))/(std::sph_bessel(0,pe*R)*std::sph_bessel(1,pe*R));
}
//经典原子核半径， L.R.B. Elton, A semi-empirical formula for the nuclear radius, Nucl. Phys. 5, 173 (1958)
double BetaDecayBranch::NuclearRadius(int A)
{
    double R = 1.121*pow(double(A),1.0/3.0)+2.426*pow(double(A),-1.0/3.0)-6.614/double(A);
    return R/197.3;
}
double BetaDecayBranch::BIG_ComplexGammaFuc(double x,double y)  //Approximation to |Gamma(x+iy)|^2  //Taken from Wilkinson, NIMA275, 378 (1989)'''
{
  double f = exp((x-1/2.)*log(x*x +y*y)-2*y*atan(y/x)-2*x+log(2*M_PI)+(1/6.)*x/(x*x+y*y));
  return f;
}
//如何区分一阶禁戒GT和一阶禁戒F???
//暂时：一切禁戒混合型均算作GT型
int BetaDecayBranch::DecayType(std::string J,std::string ParentJ,double Pi,double ParentPi,int JpiFlag,int ParentJpiFlag)
{
    
    double J1,J2,pi1,pi2,dJ,dpi;
    if(JpiFlag == 0 && ParentJpiFlag == 0)
    {
        pi1=Pi;
        pi2=ParentPi;
        std::string findstring = "/";
        int index = J.find(findstring,0);
        if(index<J.length())
        {
            J1=std::stod(J)/2.0;
            J2=std::stod(ParentJ)/2.0;
        }
        else
        {
            J1=std::stod(J);
            J2=std::stod(ParentJ);
        }
        dJ=fabs(J1-J2);
        dpi=pi1*pi2;
        if((dJ)<=1e-9 && (dpi-1)<=1e-9)
        {
            return 0;
        }
        else 
        {
            if(fabs(dJ-1)<=1e-9 && fabs(dpi-1)<=1e-9)
            {
                return 1;
            }
            else if(fabs(dJ)<=1e-9 && fabs(dpi+1)<=1e-9)
            {
                return 2;
            }
            else if(fabs(dJ-1)<=1e-9 && fabs(dpi+1)<=1e-9)
            {
                return 3;
            }
            //if((dJ-2)<=1e-9 && (dpi+1)<=1e-9)
            else
            {
                return 4;
            }
        }
    }
    //pi未知 J已知优先按Allow型
    else if((JpiFlag == 3 && (ParentJpiFlag ==0||ParentJpiFlag ==3))||(ParentJpiFlag == 3 && (JpiFlag ==0||JpiFlag ==3)))
    {
        std::string findstring = "/";
        int index = J.find(findstring,0);
        if(index<J.length())
        {
            J1=std::stod(J)/2.0;
            J2=std::stod(ParentJ)/2.0;
        }
        else
        {
            J1=std::stod(J);
            J2=std::stod(ParentJ);
        }
        dJ=fabs(J1-J2);
        if(fabs(dJ)<1e-9)
        {
            pi1=1;
            pi2=1;
            dpi=pi1*pi2;
        }
        else if(fabs(dJ-1)<1e-9)
        {
            pi1=1;
            pi2=1;
            dpi=pi1*pi2;
        }
        else
        {
            pi1=1;
            pi2=-1;
            dpi=pi1*pi2;
        }
        if((dJ)<=1e-9 && (dpi-1)<=1e-9)
        {
            return 0;
        }
        else 
        {
            if(fabs(dJ-1)<=1e-9 && fabs(dpi-1)<=1e-9)
            {
                return 1;
            }
            else if(fabs(dJ)<=1e-9 && fabs(dpi+1)<=1e-9)
            {
                return 2;
            }
            else if(fabs(dJ-1)<=1e-9 && fabs(dpi+1)<=1e-9)
            {
                return 3;
            }
            //if((dJ-2)<=1e-9 && (dpi+1)<=1e-9)
            else
            {
                return 4;
            }
        }
    }
    //Unknow与简并按GT Allow
    else
    {
        return 1;
    }
}
//1. New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
//2. Systematic Uncertainties in the Analysis of the Reactor Neutrino Anomaly
//3. Reactor antineutrino spectra and forbidden beta decays
double BetaDecayBranch::CFunc(int Z,int A,double Ee,double E0,int DecayTypeFlag,int CalTypeFlag)
{
    
    double pe=sqrt(Ee*Ee-me*me);
    double C;
    if(CalTypeFlag == 0)
    {
        
        if(DecayTypeFlag == 0) 
        {
            C=1;
        }
        if(DecayTypeFlag == 1)
        {
            C=1;
        }
        if(DecayTypeFlag == 2)
        {
            C=(pe*pe+(E0-Ee)*(E0-Ee)+2*(pe/Ee)*(pe/Ee)*Ee*(E0-Ee));
        }
        if(DecayTypeFlag == 3)
        {
            C=(pe*pe+(E0-Ee)*(E0-Ee)-(4.0/3.0)*(pe/Ee)*(pe/Ee)*Ee*(E0-Ee));
        }
        if(DecayTypeFlag == 4)
        {
            C=(pe*pe+(E0-Ee)*(E0-Ee));
        }
        if(DecayTypeFlag == 5)
        {
            C=(pe*pe+(E0-Ee)*(E0-Ee)+(2.0/3.0)*(pe/Ee)*(pe/Ee)*Ee*(E0-Ee));
        }
    }
    if(CalTypeFlag == 1)
    {
        if(DecayTypeFlag == 0) 
        {
            C=1;
        }
        if(DecayTypeFlag == 1)
        {
            C=1;
        }
        if(DecayTypeFlag == 2)
        {
            C=((E0-Ee)*(E0-Ee)+pe*pe*P12FermiLikeFunc(Z,A,Ee)+2*pe*(E0-Ee)*SP12FermiLikeFunc(Z,A,Ee));
        }
        if(DecayTypeFlag == 3)
        {
            C=((E0-Ee)*(E0-Ee)+2.0/3.0*pe*pe*P12FermiLikeFunc(Z,A,Ee)+1.0/3.0*pe*pe*P32FermiLikeFunc(Z,A,Ee)-4.0/3.0*pe*(E0-Ee)*SP12FermiLikeFunc(Z,A,Ee));
        }
        if(DecayTypeFlag == 4)
        {
            C=((E0-Ee)*(E0-Ee)+pe*pe*P32FermiLikeFunc(Z,A,Ee));
        }
        if(DecayTypeFlag == 5)
        {
            C=((E0-Ee)*(E0-Ee)+1.0/3.0*pe*pe*P12FermiLikeFunc(Z,A,Ee)+2.0/3.0*pe*pe*P32FermiLikeFunc(Z,A,Ee)+2.0/3.0*pe*(E0-Ee)*SP12FermiLikeFunc(Z,A,Ee));
        }
    }

    return C;
}
//1.Nuclear Zemach moments and finite-size corrections to allowed β decay
//2.New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
//如何区分一阶禁戒GT和一阶禁戒F???
//暂时：一切禁戒混合型均算作GT型
double BetaDecayBranch::Delta_FS(double Ee,double E0,int Z,int A,int DecayTypeFlag)
{
    double FS;
    double R = NuclearRadius(A);
    double Ev=E0-Ee;
    double r2=36.0*R/35.0;
    if(DecayTypeFlag == 0) 
    {
        FS=-1.5*Z*alpha*r2*(Ee+Ev/9+me*me/(3.0*Ee));
    }
    else
    {
        FS=-1.5*Z*alpha*r2*(Ee-Ev/27+me*me/(3.0*Ee));
    }
    return FS;
}
//Systematic Uncertainties in the Analysis of the Reactor Neutrino Anomaly
double BetaDecayBranch::Delta_WM(double Ee,double E0,int DecayTypeFlag)
{
    double pe=sqrt(Ee*Ee-me*me);
    double WM;
    if(DecayTypeFlag == 0) 
    {
        WM=0;
    }
    if(DecayTypeFlag == 1)
    {
        WM=(2.0/3.0)*((mu_v-1.0/2.0)/(mN*ga))*(Ee*(pe/Ee)*(pe/Ee)-(E0-Ee));
    }
    if(DecayTypeFlag == 2)
    {
        WM=0;
    }
    if(DecayTypeFlag == 3)
    {    
        WM=((mu_v-1.0/2.0)/(mN*ga))*((pe*pe+(E0-Ee)*(E0-Ee))*(Ee*(pe/Ee)*(pe/Ee)-(E0-Ee))+2*(pe/Ee)*(pe/Ee)*Ee*(E0-Ee)*(E0-Ee-Ee)/3.0)/(pe*pe+(E0-Ee)*(E0-Ee)-4*(pe/Ee)*(pe/Ee)*(E0-Ee)*Ee/3.0);
    }
    if(DecayTypeFlag == 4)
    {
        WM=(3.0/5.0)*((mu_v-1.0/2.0)/(mN*ga))*((pe*pe+(E0-Ee)*(E0-Ee))*(Ee*(pe/Ee)*(pe/Ee)-(E0-Ee))+2*(pe/Ee)*(pe/Ee)*Ee*(E0-Ee)*(E0-Ee-Ee)/3.0)/(pe*pe+(E0-Ee)*(E0-Ee));
    }
    if(DecayTypeFlag == 5)
    {
        WM=0;
    }
    return WM;
}
//New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
double BetaDecayBranch::e_Delta_QED(double Ee,double E0)
{
    double pe=sqrt(Ee*Ee-me*me);
    double b=pe/Ee;
    if((E0-Ee)>1e-9)
    {
        double QED=3*log(mN/me)-3.0/4.0+4*((atanh(b)/b)-1)*((E0-Ee)/(3*Ee)-3.0/2.0+log(2*(E0-Ee)/me))+4.0*LFuc(2*b/(1+b))/b+atanh(b)*(2*(1+b*b)+(E0-Ee)*(E0-Ee)/(6*Ee*Ee)-4*atanh(b))/b;   
        return alpha*QED/(2*M_PI);
    }
    else
    {
        return 0;
    }
}
//New Realization of the Conversion Calculation for Reactor Antineutrino Fluxes
double BetaDecayBranch::v_Delta_QED(double Ee,double E0)
{
    double pe=sqrt(Ee*Ee-me*me);
    double b=pe/Ee;
    double QED=3*log(mN/me)+23.0/4.0+8*LFuc(2*b/(1+b))+8*(atanh(b)/b-1)*log(2*Ee*b/me)+4*atanh(b)*((7+3*b*b)/8-2*atanh(b))/b;
    return alpha*QED/(2*M_PI);

}
double BetaDecayBranch::LFuc(double y)
{
    
    if(fabs(y)<1e-12)
    {
        return 0;
    }
    else
    {
        double value;
        TF1 *f1 = new TF1("f1","log(1-x)/x",0,y);
        value =  f1->Integral(0,y); 
        delete f1;
        return value;
    }
}
double BetaDecayBranch::BranchBetaSpec(int Z,int A,double Ee,double E0,int DecayTypeFlag,int CalTypeFlag)
{
    return sqrt(Ee*Ee-me*me)*Ee*(E0-Ee)*(E0-Ee)*FermiFunc(Z,A,Ee)*CFunc(Z,A,Ee,E0,DecayTypeFlag,CalTypeFlag)*(1+Delta_FS(Ee,E0,Z,A,DecayTypeFlag)+e_Delta_QED(Ee,E0)+Delta_WM(Ee,E0,DecayTypeFlag));
}
double BetaDecayBranch::BranchNeuSpec(int Z,int A,double Ev,double E0,int DecayTypeFlag,int CalTypeFlag)
{
    //std::cout<<E0<<" "<<FermiFunc(Z,A,(E0-Ev+me))<<" "<<CFunc(Z,A,(E0-Ev+me),E0+me,DecayTypeFlag,CalTypeFlag)<<" "<<Delta_FS((E0-Ev+me),E0+me,Z,A,DecayTypeFlag)<<" "<<v_Delta_QED((E0-Ev+me),E0+me)<<" "<<Delta_WM((E0-Ev+me),E0+me,DecayTypeFlag)<<std::endl; 
    return sqrt((E0-Ev)*(E0-Ev)-me*me)*(E0-Ev)*Ev*Ev*FermiFunc(Z,A,(E0-Ev))*CFunc(Z,A,(E0-Ev),E0,DecayTypeFlag,CalTypeFlag)*(1+Delta_FS((E0-Ev),E0,Z,A,DecayTypeFlag)+v_Delta_QED((E0-Ev),E0)+Delta_WM((E0-Ev),E0,DecayTypeFlag));

}