// Author: Qixiang Yin (yinqx@ihep.ac.cn)
// Purpose: Transfer Neu_Spectrum (E_deposited) to Neu_Spectrum (E_visible).

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"

int main()
{
    TFile* Depo_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/FY.root");
    TH1F* Depo_Spec = (TH1F*)Depo_File -> Get("h_TotalNeuSpec;1");

    TFile* Ratio_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/JUNOInputs2024_07_22.root");
    TH1D* Ratio_Graph = (TH1D*)Ratio_File -> Get("J22rc0_positronScintNL;1"); //从1.02MeV时开始记录

    double nBin = 10e2;         //沉积能谱的bin数
    double Energy_Bin = 1e-2;

    std::vector<double> E_depo; // 沉积能量
    std::vector<double> Ratio;  // 可见能量/沉积能量
    std::vector<double> E_vis;  // 可见能量

    for (int i = 1; i <= int(nBin); i++)
    {
        if (Depo_Spec -> GetBinLowEdge(i) + (1e-7) >= 1.02)
        {
            E_depo.push_back(Depo_Spec -> GetBinCenter(i));

            for (int j = 1; j <= Ratio_Graph -> GetNbinsX(); j++)
            {
                // 沉积能谱每个bin的中心能量 == Vis/Depo每个bin的左边界能量
                if (std::abs(Depo_Spec -> GetBinCenter(i) - Ratio_Graph -> GetBinLowEdge(j)) < (1e-6))
                {
                    Ratio.push_back(Ratio_Graph -> GetBinContent(j));

                    std::cout<<Depo_Spec -> GetBinCenter(i)<<' '<<Ratio_Graph -> GetBinLowEdge(j)<<std::endl;
                }
            }
        }
    }
    
    for (int i = 0; i < E_depo.size(); i++)
    {
        E_vis.push_back(Ratio[i] * E_depo[i]);
    }

    double nBin_vis = E_depo.size(); //可见能谱的bin数
    
    TH1F* h_NeuSpec_E_vis = new TH1F("h_NeuSpec_E_vis", "", int(nBin_vis), 1.02, nBin_vis * Energy_Bin);
    h_NeuSpec_E_vis -> GetXaxis() -> SetTitle("E_vis[MeV]");

    for (int i = 1; i <= int(nBin_vis); i++)
    {
        h_NeuSpec_E_vis -> SetBinContent(i,1);
    }

    TFile* outfile = new TFile("./NeuSpec_E_vis.root", "RECREATE");
    outfile -> cd();

    h_NeuSpec_E_vis -> Write();

    outfile -> Close();


    return 0;
}