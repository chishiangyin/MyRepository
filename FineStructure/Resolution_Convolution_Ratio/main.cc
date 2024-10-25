// Author: Qixiang Yin (yinqx@ihep.ac.cn)
// Purpose: Draw Ratio Spectrum between Bins in Convoluted Spectrum. 

#include <iostream>
#include <vector>
#include <string>
#include "TFile.h"
#include "TH1F.h"

int main()
{
    TFile* Convolution_Spec_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/Program/Resolution_Convolution/build/Convolution_Spec.root");
    TH1F* Convolution_Spec = (TH1F*)Convolution_Spec_File -> Get("h_Neu_Spec_E_rec;1");

    int Bin_Number = Convolution_Spec -> GetNbinsX();
    double Start = Convolution_Spec -> GetBinLowEdge(1);
    double EnergyBin = 0.01;

    TH1F* Ratio_Spec = new TH1F("Ratio_in_Convolution_Spec", "", Bin_Number, Start, Start + double(Bin_Number) * EnergyBin);
    Ratio_Spec -> GetXaxis() -> SetTitle("E[MeV]");

    for(int i = 1; i <= Bin_Number; i++)
    {
        Ratio_Spec -> SetBinContent(i,0);
    }

    for (int i = 1; i <= Bin_Number-1; i++)
    {
        if(Convolution_Spec->GetBinContent(i) != 0)
        {
            double Ratio = Convolution_Spec->GetBinContent(i+1) / Convolution_Spec->GetBinContent(i);
            Convolution_Spec -> SetBinContent(i,Ratio);
        }
    }

    TFile* outfile = new TFile("./Ratio_Convolution_Spec.root", "RECREATE");
    outfile -> cd();
    Ratio_Spec -> Write();
    outfile -> Close();
   
    return 0;
}
