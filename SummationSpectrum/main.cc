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
    
    
    if(InputType == '0' && DecayInput.size()!=0 && FissionInput.size()!=0) 
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
        }

        TFile* outfile = new TFile(OuputFile.c_str(), "RECREATE");
        outfile->cd();
        h_TotalBetaSpec->Write();
        h_TotalNeuSpec->Write();
        outfile->Close();
        return 1;
        
        

    }
    if(InputType == '1' && TestInputENSDF.size()!=0 && TestInputISO.size()!=0)
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
    }
    
    return 0;
}