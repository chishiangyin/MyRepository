#ifndef BETADECAYBRANCH_HH
#define BETADECAYBRANCH_HH
#include <iostream>
#include <cmath>
#include <string>
#include "TMath.h"
#include "TF1.h"
class BetaDecayBranch
{
    public:
        BetaDecayBranch();
        ~BetaDecayBranch();
        //费米函数
        double FermiFunc(int Z,int A,double Ee);
        double FermiFuncN(int Z,int A,double Ee,int n);
        double P32FermiLikeFunc(int Z,int A,double Ee);
        double P12FermiLikeFunc(int Z,int A,double Ee);
        double SP12FermiLikeFunc(int Z,int A,double Ee);
        //Approximation to |Gamma(x+iy)|^2  //Taken from Wilkinson, NIMA275, 378 (1989)'''
        double BIG_ComplexGammaFuc(double x,double y);
        //DecayType == 0 Allowed F ,DecayType == 1 Allowed GT, DecayType == 2 NUForbGT_0m, DecayType == 3 NUForbGT_1m, DecayType == 4 NUForbGT_2m, DecayType == 5 NUForbF_1m
        int DecayType(std::string J,std::string ParentJ,double Pi,double ParentPi,int JpiFlag,int ParentJpiFlag);
        //Shape Factor; CalTypeFlag==0: Plane Wave approximation , CalTypeFlag==1: Exact relativistic calculation;
        double CFunc(int Z,int A,double Ee,double E0,int DecayTypeFlag,int CalTypeFlag);
        //经典原子核半径
        double NuclearRadius(int A);
        //FS修正
        double Delta_FS(double Ee,double E0,int Z,int A,int DecayTypeFlag);
        //WM修正
        double Delta_WM(double Ee,double E0,int DecayTypeFlag);
        //电子QED辐射修正
        double e_Delta_QED(double Ee,double E0);
        //反中微子QED辐射修正
        double v_Delta_QED(double Ee,double E0);
        double LFuc(double y);
        //单分支Beta谱
        double BranchBetaSpec(int Z,int A,double Ee,double E0,int DecayTypeFlag,int CalTypeFlag);
        //单分支反中微子谱
        double BranchNeuSpec(int Z,int A,double Ee,double E0,int DecayTypeFlag,int CalTypeFlag);

    private:
        //电子质量
        double me;
        //fine structure constant
        double alpha;
        //中子质量
        double mN;
        //Nucleon Isovector Magnetic Moment
        double mu_v;
        //axial vector coupling constant
        double ga;
};

#endif