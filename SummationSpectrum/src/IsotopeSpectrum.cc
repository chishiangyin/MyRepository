#include "IsotopeSpectrum.hh"
#include "BetaDecayBranch.hh"
#include "Input.hh"
#include <iostream>
#include <string>
#include <numeric>
IsotopeSpectrum::IsotopeSpectrum()
{
    BetaBranchFlag = 0;
    ParentMultipleFlag = 0;
}
IsotopeSpectrum::~IsotopeSpectrum()
{
    MultParentEnergyLevel.clear();
    IsoENSDFBetaDecayData.clear();
    LevelIsoENSDFBetaDecayData.clear();
    IsoNeuSpectrum.clear();
    IsoBetaSpectrum.clear();
}
std::vector<Input::ENSDFStruct> IsotopeSpectrum::GetIsoENSDFData()
{
    return IsoENSDFBetaDecayData;
}
void IsotopeSpectrum::SetIsoENSDFData(std::string isoname,std::vector<Input::ENSDFStruct> ENSDFBetaDecayData)
{
    for(int i=0;i<ENSDFBetaDecayData.size();i++)
    {
        if(ENSDFBetaDecayData[i].isoname==isoname)
        {
            IsoENSDFBetaDecayData.push_back(ENSDFBetaDecayData[i]);
        }
    }
}

int IsotopeSpectrum::GetBetaBranchFlag()
{
    return BetaBranchFlag;
}
void IsotopeSpectrum::SetBetaBranchFlag(int a)
{
    BetaBranchFlag = a;
}
int IsotopeSpectrum::GetParentMultipleFlag()
{
    return ParentMultipleFlag;
}
void IsotopeSpectrum::SetParentMultipleFlag(int a)
{
    ParentMultipleFlag = a;
}
std::vector<double> IsotopeSpectrum::GetMultParentEnergyLevel()
{
    return MultParentEnergyLevel;
}
void IsotopeSpectrum::SetMultParentEnergyLevel(double a)
{
    MultParentEnergyLevel.push_back(a);
}
void IsotopeSpectrum::ClearIso()
{
    BetaBranchFlag = 0;
    ParentMultipleFlag = 0;
    MultParentEnergyLevel.clear();
    IsoENSDFBetaDecayData.clear();
    LevelIsoENSDFBetaDecayData.clear();
    IsoNeuSpectrum.clear();
    IsoBetaSpectrum.clear();
}
void IsotopeSpectrum::JudgeParentMultipleFlag()
{
    int flag;
    if(MultParentEnergyLevel.size()==0)
    {
        SetMultParentEnergyLevel(IsoENSDFBetaDecayData[0].ParentEnergyLevel);
    }
    for(int i=0;i<IsoENSDFBetaDecayData.size();i++)
    {
        flag = 0;
        for(int j=0;j<MultParentEnergyLevel.size();j++)
        {
            if(fabs(IsoENSDFBetaDecayData[i].ParentEnergyLevel-MultParentEnergyLevel[j])<1e-9)
            {
                flag = 1;
                break;
            }
        }
        if(flag == 0)
        {
            SetMultParentEnergyLevel(IsoENSDFBetaDecayData[i].ParentEnergyLevel);
        }
    }
    SetParentMultipleFlag(MultParentEnergyLevel.size()-1);
}
std::vector<Input::ENSDFStruct> IsotopeSpectrum::GetLevelIsoENSDFData()
{
    return  LevelIsoENSDFBetaDecayData;
}
void IsotopeSpectrum::SetLevelIsoENSDFData(int EnergyLevel)
{
    if(EnergyLevel>MultParentEnergyLevel.size())
    {
        std::cout<<"ERROR Mult Level"<<std::endl;
        return;
    }
    std::sort(MultParentEnergyLevel.begin(),MultParentEnergyLevel.end());
    double energy = MultParentEnergyLevel[EnergyLevel];
    for(int i=0;i<IsoENSDFBetaDecayData.size();i++)
    {
        if(fabs(energy-IsoENSDFBetaDecayData[i].ParentEnergyLevel)<1e-9)
        {
            LevelIsoENSDFBetaDecayData.push_back(IsoENSDFBetaDecayData[i]);
        }
    }
}   
void IsotopeSpectrum::JudgeBetaBranchFlag()
{
    double Beta_Branch_sum=0;
    for(int i=0;i<LevelIsoENSDFBetaDecayData.size();i++)
    {
        Beta_Branch_sum+=LevelIsoENSDFBetaDecayData[i].BranchStrength;
    }
    if(Beta_Branch_sum>0.99 && Beta_Branch_sum<1.01)
    {
        SetBetaBranchFlag(0);
    }
    else
    {
        SetBetaBranchFlag(1);
    }

}
void IsotopeSpectrum::BranchNormal()
{
    double Beta_Branch_sum=0;
    for(int i=0;i<LevelIsoENSDFBetaDecayData.size();i++)
    {
        Beta_Branch_sum+=LevelIsoENSDFBetaDecayData[i].BranchStrength;
    }
    if(Beta_Branch_sum<1e-5)
    {
        return;
    }
    for(int i=0;i<LevelIsoENSDFBetaDecayData.size();i++)
    {
        LevelIsoENSDFBetaDecayData[i].BranchStrength/=Beta_Branch_sum;
    }

}
void IsotopeSpectrum::SetIsoBetaSpec(int CalTypeFlag,int BranchNormalFlag,double EnergyBin)
{
    BetaDecayBranch BetaDecay;
    if(BranchNormalFlag==1)
    {
        BranchNormal();
    }
    std::vector<double> BranchBetaSpectrum;
    for(int i=0;i<LevelIsoENSDFBetaDecayData.size();i++)
    {
        double me=0.511099895;
        BranchBetaSpectrum.clear();
        int DT = BetaDecay.DecayType(LevelIsoENSDFBetaDecayData[i].J,LevelIsoENSDFBetaDecayData[i].ParentJ,LevelIsoENSDFBetaDecayData[i].Pi,LevelIsoENSDFBetaDecayData[i].ParentPi,LevelIsoENSDFBetaDecayData[i].JpiFlag,LevelIsoENSDFBetaDecayData[i].ParentJpiFlag);
        double E0 = (LevelIsoENSDFBetaDecayData[i].Qvalue+LevelIsoENSDFBetaDecayData[i].ParentEnergyLevel-LevelIsoENSDFBetaDecayData[i].EnergyLevel);
        E0+=me;
        for(int j=0;j<int((LevelIsoENSDFBetaDecayData[i].Qvalue+LevelIsoENSDFBetaDecayData[i].ParentEnergyLevel-LevelIsoENSDFBetaDecayData[i].EnergyLevel)/EnergyBin)+1;j++)
        {
            double Te=double(j+0.5)*EnergyBin;      
            double Ee=Te+me;
            if(Ee>1e-9 && Ee<E0-1e-9)
            {
                BranchBetaSpectrum.push_back(BetaDecay.BranchBetaSpec(LevelIsoENSDFBetaDecayData[i].Z,LevelIsoENSDFBetaDecayData[i].A,Ee,E0,DT,CalTypeFlag));
            }
            else
            {
                BranchBetaSpectrum.push_back(0);
            }

        }
        double norm = std::accumulate(BranchBetaSpectrum.begin(), BranchBetaSpectrum.end(), 0.0);
        for(int j=0;j<BranchBetaSpectrum.size();j++)
        {
            BranchBetaSpectrum[j]/= (norm * EnergyBin);
            if(j>=IsoBetaSpectrum.size())
            {
                IsoBetaSpectrum.push_back(BranchBetaSpectrum[j]*LevelIsoENSDFBetaDecayData[i].BranchStrength);
            }
            else
            {
                IsoBetaSpectrum[j]+=BranchBetaSpectrum[j]*LevelIsoENSDFBetaDecayData[i].BranchStrength;
            }
        }
    }


}
void IsotopeSpectrum::SetIsoNeuSpec(int CalTypeFlag,int BranchNormalFlag,double EnergyBin)
{
    BetaDecayBranch BetaDecay;
    if(BranchNormalFlag==1)
    {
        BranchNormal();
    }
    std::vector<double> BranchNeuSpectrum;
    //std::cout<<"LEVEL SIZE:"<<LevelIsoENSDFBetaDecayData.size()<<std::endl;
    for(int i=0;i<LevelIsoENSDFBetaDecayData.size();i++)
    {
        BranchNeuSpectrum.clear();
        int DT = BetaDecay.DecayType(LevelIsoENSDFBetaDecayData[i].J,LevelIsoENSDFBetaDecayData[i].ParentJ,LevelIsoENSDFBetaDecayData[i].Pi,LevelIsoENSDFBetaDecayData[i].ParentPi,LevelIsoENSDFBetaDecayData[i].JpiFlag,LevelIsoENSDFBetaDecayData[i].ParentJpiFlag);
        double E0 = (LevelIsoENSDFBetaDecayData[i].Qvalue+LevelIsoENSDFBetaDecayData[i].ParentEnergyLevel-LevelIsoENSDFBetaDecayData[i].EnergyLevel);
        double me=0.511099895;
        E0+=me;
        for(int j=0;j<int((LevelIsoENSDFBetaDecayData[i].Qvalue+LevelIsoENSDFBetaDecayData[i].ParentEnergyLevel-LevelIsoENSDFBetaDecayData[i].EnergyLevel)/EnergyBin)+1;j++)
        {
            double Ev=double(j+0.5)*EnergyBin;
            if(Ev>1e-9 && Ev<(E0-me)-1e-9)
            {
                BranchNeuSpectrum.push_back(BetaDecay.BranchNeuSpec(LevelIsoENSDFBetaDecayData[i].Z,LevelIsoENSDFBetaDecayData[i].A,Ev,E0,DT,CalTypeFlag));
                //std::cout<<LevelIsoENSDFBetaDecayData[i].BranchStrength<<" "<<LevelIsoENSDFBetaDecayData[i].EnergyLevel<<" "<<BetaDecay.BranchNeuSpec(LevelIsoENSDFBetaDecayData[i].Z,LevelIsoENSDFBetaDecayData[i].A,Ev,E0,DT,CalTypeFlag)<<std::endl;
            }
            else
            {
                BranchNeuSpectrum.push_back(0);
            }

        }
        double norm = std::accumulate(BranchNeuSpectrum.begin(), BranchNeuSpectrum.end(), 0.0);
        if(norm==0)
        {
            std::cout<<"ERROR norm"<<std::endl;
        }
        for(int j=0;j<BranchNeuSpectrum.size();j++)
        {
            BranchNeuSpectrum[j]/= (norm * EnergyBin);
            if(j>=IsoNeuSpectrum.size())
            {
                IsoNeuSpectrum.push_back(BranchNeuSpectrum[j]*LevelIsoENSDFBetaDecayData[i].BranchStrength);
            }
            else
            {
                IsoNeuSpectrum[j]+=BranchNeuSpectrum[j]*LevelIsoENSDFBetaDecayData[i].BranchStrength;
            }
        }
    }
}
std::vector<double> IsotopeSpectrum::GetIsoBetaSpec()
{
    return IsoBetaSpectrum;
}
std::vector<double> IsotopeSpectrum::GetIsoNeuSpec()
{
    return IsoNeuSpectrum;
}
int IsotopeSpectrum::GetNIso()
{
    return IsoENSDFBetaDecayData.size();
}
int IsotopeSpectrum::GetLevelNIso()
{
    return LevelIsoENSDFBetaDecayData.size();
}