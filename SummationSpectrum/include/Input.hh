#ifndef INPUT_HH
#define INPUT_HH
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <sstream>
#include <fstream>
class Input
{
    public:
        Input();
        ~Input();
        //Flag:0:Jpi唯一  1：Jpi未知  2：Jpi不唯一    3：J确定 pi 不确定
        //J：只要有信息就记录，唯一记录数值，不唯一全信息记录，包括pi的信息
        //pi: 1,-1:为pi确定时的数据记录   2：J确定pi不确定   3：说明有兼并，pi信息由J记录  4:说明Jpi未知
        struct ENSDFStruct 
        {
            std::string isoname;
            int A;
            int Z;
            double Qvalue;
            double BranchStrength;
            double EnergyLevel;
            std::string J;
            double Pi;
            double ParentEnergyLevel;
            std::string ParentJ;
            double ParentPi;
            int JpiFlag;
            int ParentJpiFlag;
        };
        struct ISOStruct
        {
            std::string isoname;
            double activity;
        };
        struct ISOFYStruct
        {
            std::string Parentisoname;
            int ParentA;
            int ParentZ;
            std::string isoname;
            int A;
            int Z;
            int NEnergyLevel;
            double thermal_fy;
            double thermal_fy_unc;
            double fast_fy;
            double fast_fy_unc;
            double fourteen_mev_fy;
            double fourteen_mev_fy_unc;
        };
        std::vector<Input::ISOFYStruct> GetFYList();
        void SetFYList(Input::ISOFYStruct iso);
        std::vector<Input::ISOStruct> GetIsoList();
        void SetIsoList(Input::ISOStruct iso);
        std::vector<Input::ENSDFStruct> GetENSDFData();
        void SetENSDFData(Input::ENSDFStruct ENSDF);
        void SetIsoFileName(std::string filename);
        void SetENSDFFileName(std::string filename);
        void SetFYFileName(std::string filename);
        void ReadENSDFInputFile();
        void ReadFYInputFile();
        void ReadIsoInputFile();
    private:
        std::vector<Input::ISOFYStruct> IsoFYList;
        std::vector<Input::ISOStruct> IsoActivityList;
        std::vector<Input::ENSDFStruct> ENSDFBetaDecayData;
        std::string isofilename;
        std::string ENSDFfilename;
        std::string FYfilename;
};
#endif