// Modifying by Qixiang Yin (yinqx@ihep.ac.cn)
// New: 可以产生所有裂变分支的产额谱以及Ratio谱，输出在8MeV处导致精细结构的核素（Version 2）

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "BetaDecayBranch.hh"
#include "BetaDecayInput.hh"
#include "FissionInput.hh"
#include "InitializeInput.hh"
#include "IsotopeSpectrum.hh"
#include "Input.hh"
#include <string>
#include <algorithm>

int main(int argc, char* argv[]) 
{
    //ShapeFactorCalFlag==0:Plane Wave approximation, ShapeFactorCalFlag==1: Exact relativistic calculation;
    std::string ShapeFactorCalFlag;
    //InputType==0:Use Fission Input, InputType==1:Use Isotope Input;
    std::string InputType;
    //NeutronType==0:Thermal, NeutronType==1:Fast, NeutronType==2:14MeV;
    std::string NeutronType;
    std::string IsotopeInput;
    std::string FissionInput;
    std::string DecayInput;
    std::string OuputFile;
    std::string InitInput;
    std::string TestInputENSDF;
    std::string TestInputISO;
    std::string EnergyBin;
    std::string nBin;
    //Normalized beta decay branching ratio of isotope, 0 No, 1 Yes 
    std::string BranchNormalFlag;
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-ShapeFactorCalFlag") {
            if (i + 1 < argc) {
                ShapeFactorCalFlag = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-InputType") {
            if (i + 1 < argc) {
                InputType = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-IsotopeInput") {
            if (i + 1 < argc) {
                IsotopeInput = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-FissionInput") {
            if (i + 1 < argc) {
                FissionInput = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-NeutronType") {
            if (i + 1 < argc) {
                NeutronType = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-DecayInput") {
            if (i + 1 < argc) {
                DecayInput = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-OuputFile") {
            if (i + 1 < argc) {
                OuputFile = argv[i + 1];
            }
        }
        
        if (std::string(argv[i]) == "-TestInputENSDF") {
            if (i + 1 < argc) {
                TestInputENSDF = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-TestInputISO") {
            if (i + 1 < argc) {
                TestInputISO = argv[i + 1];
            }
        }
        //Unit:MeV
        if (std::string(argv[i]) == "-EnergyBin") {
            if (i + 1 < argc) {
                EnergyBin = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-nBin") {
            if (i + 1 < argc) {
                nBin = argv[i + 1];
            }
        }
        if (std::string(argv[i]) == "-BranchNormalFlag") {
            if (i + 1 < argc) {
                BranchNormalFlag = argv[i + 1];
            }
        }
    }
    if(OuputFile.size()==0)
    {
        OuputFile="output.root";
    }
    if(ShapeFactorCalFlag.size()==0)
    {
        ShapeFactorCalFlag="0";
    }
    if(BranchNormalFlag.size()==0)
    {
        BranchNormalFlag="1";
    }
    if(EnergyBin.size()==0)
    {
        EnergyBin="1e-2";
    }
    if(nBin.size()==0)
    {
        nBin="10e2";
    }
    
    //Select FineStructure Isotope
    struct FineStructure
    {
        std::string Iso_Name;
        double Neu_Spec_809 = 0.0;
        double Neu_Spec_810 = 0.0;
        double Ratio = Neu_Spec_810/Neu_Spec_809;
    };
    
    if(InputType == '0' && DecayInput.size()!=0 && FissionInput.size()!=0) 
    {
        
        TH1F* h_TotalBetaSpec = new TH1F("h_TotalBetaSpec","",int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
        h_TotalBetaSpec->GetXaxis()->SetTitle("E[MeV]"); 
        TH1F* h_TotalNeuSpec = new TH1F("h_TotalNeuSpec","",int(stod(nBin)),0,stod(nBin)*stof(EnergyBin));
        h_TotalNeuSpec->GetXaxis()->SetTitle("E[MeV]");

        //Ratio
        TH1F* h_TotalBetaRatioSpec = new TH1F("h_TotalBetaRatioSpec","",int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
        h_TotalBetaRatioSpec->GetXaxis()->SetTitle("E[MeV]");
        TH1F* h_TotalNeuRatioSpec = new TH1F("h_TotalNeuRatioSpec","",int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
        h_TotalNeuRatioSpec->GetXaxis()->SetTitle("E[MeV]");

        for(int i=0;i<int(stod(nBin));i++)
        {
            h_TotalBetaSpec->SetBinContent(i+1,0);
            h_TotalNeuSpec->SetBinContent(i+1,0);
        }

        //Ratio
        for(int i=0; i<int(stod(nBin)); i++)
        {
            h_TotalBetaRatioSpec->SetBinContent(i+1,0);
            h_TotalNeuRatioSpec->SetBinContent(i+1,0);
        } 

        std::cout<<OuputFile<<std::endl;
        std::cout<<ShapeFactorCalFlag<<std::endl;
        std::cout<<FissionInput<<std::endl;
        std::cout<<DecayInput<<std::endl;
        std::cout<<BranchNormalFlag<<std::endl;
        std::cout<<EnergyBin<<std::endl;
        ///////////////////////////////
        Input Input1;
        Input1.SetENSDFFileName(DecayInput);
        Input1.SetFYFileName(FissionInput);
        Input1.ReadENSDFInputFile();
        Input1.ReadFYInputFile();
        std::vector<Input::ENSDFStruct> ENSDFList;
        std::vector<Input::ISOFYStruct> FYList;
        ENSDFList = Input1.GetENSDFData();
        FYList = Input1.GetFYList();

        //Select FineStructure Isotope
        std::ofstream output("Selected_FineStructure.txt");
        std::vector<FineStructure> FineStructure_8MeV;

        for(int i=0;i<FYList.size();i++)
        {
        
            int energylevel = FYList[i].NEnergyLevel;
            std::string sisoname;
            sisoname = FYList[i].isoname;
            std::cout<<"Isotope Name:"<<sisoname<<std::endl;
            IsotopeSpectrum IsotopeSpectrum1;
            IsotopeSpectrum1.ClearIso();
            IsotopeSpectrum1.SetIsoENSDFData(sisoname,ENSDFList);
            if(IsotopeSpectrum1.GetNIso()==0)
            {
                //std::cout<<"ERROR No Isotope in ENSDF"<<std::endl;
                std::cout<<"ERROR No Beta Isotope:"<<sisoname<<"  in ENSDF"<<std::endl;
                continue;
            }
            IsotopeSpectrum1.JudgeParentMultipleFlag();
            //需要进一步修正<0的能级？
            if(energylevel<=0)
            {
                IsotopeSpectrum1.SetLevelIsoENSDFData(0);
            }
            else if(IsotopeSpectrum1.GetParentMultipleFlag()<energylevel)
            {
                IsotopeSpectrum1.SetLevelIsoENSDFData(IsotopeSpectrum1.GetParentMultipleFlag());
            }
            else
            {
                IsotopeSpectrum1.SetLevelIsoENSDFData(energylevel);
            }
            if(IsotopeSpectrum1.GetLevelNIso()==0)
            {
                std::cout<<"WANRING No Level Isotope in ENSDF"<<std::endl;
                continue;
            }

            IsotopeSpectrum1.JudgeBetaBranchFlag();
            IsotopeSpectrum1.SetIsoBetaSpec(stoi(ShapeFactorCalFlag),stoi(BranchNormalFlag),stod(EnergyBin));
            IsotopeSpectrum1.SetIsoNeuSpec(stoi(ShapeFactorCalFlag),stoi(BranchNormalFlag),stod(EnergyBin));
            std::vector<double> IsoBetaSpec;
            std::vector<double> IsoNeuSpec;
            IsoBetaSpec = IsotopeSpectrum1.GetIsoBetaSpec();
            IsoNeuSpec = IsotopeSpectrum1.GetIsoNeuSpec();

            //Selected FineStructure Isotope
            double Neu_Spec_809 = 0.0;
            double Neu_Spec_810 = 0.0;

            for(int j=0;j<int(stod(nBin));j++)
            {
                
                if(NeutronType=='0')
                {

                
                    if(j<IsoBetaSpec.size())
                    {
                        h_TotalBetaSpec->SetBinContent(j+1,h_TotalBetaSpec->GetBinContent(j+1)+IsoBetaSpec[j]*FYList[i].thermal_fy);
                    }
                    if(j<IsoNeuSpec.size())
                    {
                        h_TotalNeuSpec->SetBinContent(j+1,h_TotalNeuSpec->GetBinContent(j+1)+IsoNeuSpec[j]*FYList[i].thermal_fy);

                        //Selected FineStructure Isotope
                        if(j+1 == 809)
                        {
                            Neu_Spec_809 = IsoNeuSpec[j]*FYList[i].thermal_fy;
                        }
                        if(j+1 == 810)
                        {
                            Neu_Spec_810 = IsoNeuSpec[j]*FYList[i].thermal_fy;
                        }
                    }
                }

                if(NeutronType=='1')
                {

                
                    if(j<IsoBetaSpec.size())
                    {
                        h_TotalBetaSpec->SetBinContent(j+1,h_TotalBetaSpec->GetBinContent(j+1)+IsoBetaSpec[j]*FYList[i].fast_fy);
                    }
                    if(j<IsoNeuSpec.size())
                    {
                        h_TotalNeuSpec->SetBinContent(j+1,h_TotalNeuSpec->GetBinContent(j+1)+IsoNeuSpec[j]*FYList[i].fast_fy);
                    }

                }

                if(NeutronType=='2')
                {

                
                    if(j<IsoBetaSpec.size())
                    {
                        h_TotalBetaSpec->SetBinContent(j+1,h_TotalBetaSpec->GetBinContent(j+1)+IsoBetaSpec[j]*FYList[i].fourteen_mev_fy);
                    }
                    if(j<IsoNeuSpec.size())
                    {
                        h_TotalNeuSpec->SetBinContent(j+1,h_TotalNeuSpec->GetBinContent(j+1)+IsoNeuSpec[j]*FYList[i].fourteen_mev_fy);
                    }

                }
                
            }

            if(Neu_Spec_809 != 0.0)
            {
                FineStructure_8MeV.push_back({sisoname,Neu_Spec_809,Neu_Spec_810});
            }
        }

        //Selected FineStructure Isotope
        for(int i=0; i<FineStructure_8MeV.size(); i++)
        {
            for(int j=0; j<FineStructure_8MeV.size()-i-1; j++)
            {
                if(FineStructure_8MeV[j].Ratio > FineStructure_8MeV[j+1].Ratio)
                {
                    std::swap(FineStructure_8MeV[j],FineStructure_8MeV[j+1]);
                }
            }
        }

        output << "----IsoName----809----810----Ratio----" << std::endl;

        for(int i=0; i<FineStructure_8MeV.size(); i++)
        {
            output << FineStructure_8MeV[i].Iso_Name << ' ' << FineStructure_8MeV[i].Neu_Spec_809 << ' ' << FineStructure_8MeV[i].Neu_Spec_810 << ' ' << FineStructure_8MeV[i].Ratio << std::endl;
        }

        double Neu_Spec_809_Sum = 0.0;
        double Neu_Spec_810_Sum = 0.0;

        for(int i=0; i<FineStructure_8MeV.size(); i++)
        {
            Neu_Spec_809_Sum += FineStructure_8MeV[i].Neu_Spec_809;
        }

        for(int i=0; i<FineStructure_8MeV.size(); i++)
        {
            Neu_Spec_810_Sum += FineStructure_8MeV[i].Neu_Spec_810;
        }

        output << "----------Total Ratio-----------" << std::endl;

        output << "----809----810----Ratio----" << std::endl;

        output << Neu_Spec_809_Sum << ' ' << Neu_Spec_810_Sum << ' ' << Neu_Spec_810_Sum/Neu_Spec_809_Sum << std::endl;

        //Ratio
        for(int i=0; i<int(stod(nBin)); i++)
        {
            if(NeutronType == '0')
            {
                h_TotalBetaRatioSpec->SetBinContent(i+1,h_TotalBetaSpec->GetBinContent(i+2)/h_TotalBetaSpec->GetBinContent(i+1));
                h_TotalNeuRatioSpec->SetBinContent(i+1,h_TotalNeuSpec->GetBinContent(i+2)/h_TotalNeuSpec->GetBinContent(i+1));
            }
        }

        //Branch FY Spectrum
        int Flag=0;
        for(int i=0; i<FYList.size(); i++)
        {
            int energylevel = FYList[i].NEnergyLevel;
            std::string sisoname;
            sisoname = FYList[i].isoname;
            IsotopeSpectrum IsotopeSpectrum1;
            IsotopeSpectrum1.ClearIso();
            IsotopeSpectrum1.SetIsoENSDFData(sisoname,ENSDFList);
            if(IsotopeSpectrum1.GetNIso()==0)
            {
                continue;
            }
            else
            {
                Flag++;
            }
        }

        TH1F* h_BetaSpec[Flag];
        TH1F* h_NeuSpec[Flag];

        //Ratio
        TH1F* h_BetaSpecRatio[Flag];
        TH1F* h_NeuSpecRatio[Flag];

        Flag=0;
        for(int i=0; i<FYList.size(); i++)
        {
            int energylevel = FYList[i].NEnergyLevel;
            std::string sisoname;
            sisoname = FYList[i].isoname;
            IsotopeSpectrum IsotopeSpectrum1;
            IsotopeSpectrum1.ClearIso();
            IsotopeSpectrum1.SetIsoENSDFData(sisoname,ENSDFList);
            if(IsotopeSpectrum1.GetNIso()==0)
            {
                continue;
            }
            std::string betaname;
            std::string neuname;
            std::string EnergyLevel_str = std::to_string(energylevel);
            betaname = "h_BetaSpec_";
            neuname = "h_NeuSpec_";
            betaname+=sisoname;
            betaname+='_';
            betaname+=EnergyLevel_str;
            neuname+=sisoname;
            neuname+='_';
            neuname+=EnergyLevel_str;
            h_BetaSpec[Flag] = new TH1F(betaname.c_str(),betaname.c_str(),int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
            h_NeuSpec[Flag] = new TH1F(neuname.c_str(),neuname.c_str(),int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
            h_BetaSpec[Flag]->GetXaxis()->SetTitle("E[MeV]");
            h_NeuSpec[Flag]->GetXaxis()->SetTitle("E[MeV]");

            //Ratio
            std::string betaname_2 = betaname+="_Ratio";
            std::string neuname_2 = neuname+="_Ratio";
            h_BetaSpecRatio[Flag] = new TH1F(betaname_2.c_str(),betaname_2.c_str(),int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
            h_NeuSpecRatio[Flag] = new TH1F(neuname_2.c_str(),neuname_2.c_str(),int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
            h_BetaSpecRatio[Flag]->GetXaxis()->SetTitle("E[MeV]");
            h_NeuSpecRatio[Flag]->GetXaxis()->SetTitle("E[MeV]");

            std::cout<<"Isotope:"<<sisoname<<' '<<"EnergyLevel:"<<EnergyLevel_str<<std::endl;
            for(int j=0; j<int(stod(nBin)); j++)
            {
                h_BetaSpec[Flag]->SetBinContent(j+1,0);
                h_NeuSpec[Flag]->SetBinContent(j+1,0);

                //Ratio
                h_BetaSpecRatio[Flag]->SetBinContent(j+1,0);
                h_NeuSpecRatio[Flag]->SetBinContent(j+1,0);
            }
            IsotopeSpectrum1.JudgeParentMultipleFlag();
            if(energylevel<=0)
            {
                IsotopeSpectrum1.SetLevelIsoENSDFData(0);
            }
            else if(IsotopeSpectrum1.GetParentMultipleFlag()<energylevel)
            {
                IsotopeSpectrum1.SetLevelIsoENSDFData(IsotopeSpectrum1.GetParentMultipleFlag());
            }
            else
            {
                IsotopeSpectrum1.SetLevelIsoENSDFData(energylevel);
            }
            IsotopeSpectrum1.JudgeBetaBranchFlag();
            IsotopeSpectrum1.SetIsoBetaSpec(stoi(ShapeFactorCalFlag),stoi(BranchNormalFlag),stod(EnergyBin));
            IsotopeSpectrum1.SetIsoNeuSpec(stoi(ShapeFactorCalFlag),stoi(BranchNormalFlag),stod(EnergyBin));
            std::vector<double> IsoBetaSpec;
            std::vector<double> IsoNeuSpec;
            IsoBetaSpec = IsotopeSpectrum1.GetIsoBetaSpec();
            IsoNeuSpec = IsotopeSpectrum1.GetIsoNeuSpec();
            for(int j=0;j<int(stod(nBin));j++)
            {
                if(NeutronType=='0')
                {
                    if(j<IsoBetaSpec.size())
                    {
                        h_BetaSpec[Flag]->SetBinContent(j+1,h_BetaSpec[Flag]->GetBinContent(j+1)+IsoBetaSpec[j]*FYList[i].thermal_fy);
                    }
                    if(j<IsoNeuSpec.size())
                    {
                        h_NeuSpec[Flag]->SetBinContent(j+1,h_NeuSpec[Flag]->GetBinContent(j+1)+IsoNeuSpec[j]*FYList[i].thermal_fy);
                    }

                }
            }
            Flag++; 
        }

        //Ratio
        for(int i=0; i<Flag; i++)
        {
            for(int j=0; j<int(stod(nBin)); j++)
            {
                if(h_BetaSpec[i]->GetBinContent(j+1) != 0)
                {
                    h_BetaSpecRatio[i]->SetBinContent(j+1,h_BetaSpec[i]->GetBinContent(j+2)/h_BetaSpec[i]->GetBinContent(j+1));
                }
                if(h_NeuSpec[i]->GetBinContent(j+1) != 0)
                {
                    h_NeuSpecRatio[i]->SetBinContent(j+1,h_NeuSpec[i]->GetBinContent(j+2)/h_NeuSpec[i]->GetBinContent(j+1));
                }
                if(h_BetaSpec[i]->GetBinContent(j+1) == 0)
                {
                    h_BetaSpecRatio[i]->SetBinContent(j+1,0);
                }
                if(h_NeuSpec[i]->GetBinContent(j+1) == 0)
                {
                    h_NeuSpecRatio[i]->SetBinContent(j+1,0);
                }
            }
        }

        TFile* outfile = new TFile(OuputFile.c_str(), "RECREATE");
        outfile->cd();
        h_TotalBetaSpec->Write();
        h_TotalNeuSpec->Write();
        
        //ratio
        h_TotalBetaRatioSpec->Write();
        h_TotalNeuRatioSpec->Write();

        //branch FY spectrum
        for(int i=0; i<Flag; i++)
        {
            h_BetaSpec[i]->Write();
            h_NeuSpec[i]->Write();

            //ratio
            h_BetaSpecRatio[i]->Write();
            h_NeuSpecRatio[i]->Write();
        }

        outfile->Close();
        return 1;
        
        

    }
    /* if(InputType == '1' && TestInputENSDF.size()!=0 && TestInputISO.size()!=0)
    {
        TH1F* h_TotalBetaSpec = new TH1F("h_TotalBetaSpec","",int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
        h_TotalBetaSpec->GetXaxis()->SetTitle("E[MeV]"); 
        TH1F* h_TotalNeuSpec = new TH1F("h_TotalNeuSpec","",int(stod(nBin)),0,stod(nBin)*stof(EnergyBin));
        h_TotalNeuSpec->GetXaxis()->SetTitle("E[MeV]");
        for(int i=0;i<int(stod(nBin));i++)
        {
            h_TotalBetaSpec->SetBinContent(i+1,0);
            h_TotalNeuSpec->SetBinContent(i+1,0);
        }
        std::cout<<OuputFile<<std::endl;
        std::cout<<ShapeFactorCalFlag<<std::endl;
        std::cout<<TestInputENSDF<<std::endl;
        std::cout<<TestInputISO<<std::endl;
        std::cout<<BranchNormalFlag<<std::endl;
        std::cout<<EnergyBin<<std::endl;
        Input Input1;
        Input1.SetENSDFFileName(TestInputENSDF);
        Input1.SetIsoFileName(TestInputISO);
        Input1.ReadENSDFInputFile();
        Input1.ReadIsoInputFile();
        std::vector<Input::ENSDFStruct> ENSDFList;
        std::vector<Input::ISOStruct> ISOList;
        ENSDFList = Input1.GetENSDFData();
        ISOList = Input1.GetIsoList();
        TH1F* h_BetaSpec[ISOList.size()];
        TH1F* h_NeuSpec[ISOList.size()];
        for(int i=0;i<ISOList.size();i++)
        {
            std::string betaname;
            std::string neuname;
            betaname = "h_BetaSpec_";
            neuname = "h_NeuSpec_";
            betaname+=ISOList[i].isoname;
            neuname+=ISOList[i].isoname;
            h_BetaSpec[i] = new TH1F(betaname.c_str(),betaname.c_str(),int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
            h_NeuSpec[i] = new TH1F(neuname.c_str(),neuname.c_str(),int(stod(nBin)),0,stod(nBin)*stod(EnergyBin));
            h_BetaSpec[i]->GetXaxis()->SetTitle("E[MeV]");
            h_NeuSpec[i]->GetXaxis()->SetTitle("E[MeV]");
            for(int j=0;j<int(stod(nBin));j++)
            {
                h_BetaSpec[i]->SetBinContent(j+1,0);
                h_NeuSpec[i]->SetBinContent(j+1,0);
            } 
        }

        for(int i=0;i<ISOList.size();i++)
        {
            //need update
            int energylevel = 0;
            std::string sisoname;
            if(ISOList[i].isoname.at(0)=='m')
            {
                energylevel = 1;
                sisoname = ISOList[i].isoname.substr(1);
            }
            else
            {
                sisoname = ISOList[i].isoname;
            }
            std::cout<<"Isotope Name:"<<ISOList[i].isoname<<std::endl;
            std::cout<<"Isotope Name:"<<sisoname<<std::endl;
            IsotopeSpectrum IsotopeSpectrum1;
            IsotopeSpectrum1.ClearIso();
            IsotopeSpectrum1.SetIsoENSDFData(sisoname,ENSDFList);
            if(IsotopeSpectrum1.GetNIso()==0)
            {
                //std::cout<<"ERROR No Isotope in ENSDF"<<std::endl;
                std::cout<<"ERROR No Isotope:"<<sisoname<<"  in ENSDF"<<std::endl;
                continue;
            }
            IsotopeSpectrum1.JudgeParentMultipleFlag();
            //need update
            IsotopeSpectrum1.SetLevelIsoENSDFData(energylevel);
            if(IsotopeSpectrum1.GetLevelNIso()==0)
            {
                std::cout<<"WANRING No Level Isotope in ENSDF"<<std::endl;
                continue;
            }
            
            IsotopeSpectrum1.JudgeBetaBranchFlag();
            IsotopeSpectrum1.SetIsoBetaSpec(stoi(ShapeFactorCalFlag),stoi(BranchNormalFlag),stod(EnergyBin));
            IsotopeSpectrum1.SetIsoNeuSpec(stoi(ShapeFactorCalFlag),stoi(BranchNormalFlag),stod(EnergyBin));
            std::vector<double> IsoBetaSpec;
            std::vector<double> IsoNeuSpec;
            IsoBetaSpec = IsotopeSpectrum1.GetIsoBetaSpec();
            IsoNeuSpec = IsotopeSpectrum1.GetIsoNeuSpec();
            for(int j=0;j<int(stod(nBin));j++)
            {
                if(j<IsoBetaSpec.size())
                {
                    h_TotalBetaSpec->SetBinContent(j+1,h_TotalBetaSpec->GetBinContent(j+1)+IsoBetaSpec[j]*ISOList[i].activity);
                }
                if(j<IsoNeuSpec.size())
                {
                    h_TotalNeuSpec->SetBinContent(j+1,h_TotalNeuSpec->GetBinContent(j+1)+IsoNeuSpec[j]*ISOList[i].activity);
                }
                
            }
            for(int j=0;j<int(stod(nBin));j++)
            {
                if(j<IsoBetaSpec.size())
                {
                    h_BetaSpec[i]->SetBinContent(j+1,h_BetaSpec[i]->GetBinContent(j+1)+IsoBetaSpec[j]*ISOList[i].activity);
                }
                if(j<IsoNeuSpec.size())
                {
                    h_NeuSpec[i]->SetBinContent(j+1,h_NeuSpec[i]->GetBinContent(j+1)+IsoNeuSpec[j]*ISOList[i].activity);
                }   
            }
        }
        TFile* outfile = new TFile(OuputFile.c_str(), "RECREATE");
        outfile->cd();
        h_TotalBetaSpec->Write();
        h_TotalNeuSpec->Write();
        for(int i=0;i<ISOList.size();i++)
        {
            h_BetaSpec[i]->Write();
            h_NeuSpec[i]->Write();
        }
        outfile->Close();
        return 1;
    } */
    
    return 0;
}