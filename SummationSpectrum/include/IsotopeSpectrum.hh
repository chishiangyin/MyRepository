#ifndef ISOTOPESPECTRUM_HH
#define ISOTOPESPECTRUM_HH
#include "Input.hh"
#include <iostream>
#include <string>
#include <vector>
class IsotopeSpectrum
{
    public:
        IsotopeSpectrum();
        ~IsotopeSpectrum();
        std::vector<Input::ENSDFStruct> GetIsoENSDFData();
        void SetIsoENSDFData(std::string isoname,std::vector<Input::ENSDFStruct> ENSDFBetaDecayData);
        int GetBetaBranchFlag();
        void SetBetaBranchFlag(int a);
        int GetParentMultipleFlag();
        void SetParentMultipleFlag(int a);
        std::vector<double> GetMultParentEnergyLevel();
        void SetMultParentEnergyLevel(double a);
        void JudgeParentMultipleFlag();
        std::vector<Input::ENSDFStruct> GetLevelIsoENSDFData();
        void SetLevelIsoENSDFData(int EnergyLevel);
        void JudgeBetaBranchFlag();
        void BranchNormal();
        void ClearIso();
        void SetIsoBetaSpec(int CalTypeFlag,int BranchNormalFlag,double EnergyBin);
        void SetIsoNeuSpec(int CalTypeFlag,int BranchNormalFlag,double EnergyBin);
        std::vector<double> GetIsoBetaSpec();
        std::vector<double> GetIsoNeuSpec();
        int GetNIso();
        int GetLevelNIso();
        
    private:
        std::vector<Input::ENSDFStruct> IsoENSDFBetaDecayData;
        std::vector<Input::ENSDFStruct> LevelIsoENSDFBetaDecayData;
        //Branch归一化:0归一化，1非归一化
        int BetaBranchFlag;
        //母核多能级:0无，>0有
        int ParentMultipleFlag;
        //母核能级
        std::vector<double> MultParentEnergyLevel;
        //Branch Spec
        // 
        //Iso Spec
        std::vector<double> IsoNeuSpectrum;
        std::vector<double> IsoBetaSpectrum;
    
        
};  
#endif