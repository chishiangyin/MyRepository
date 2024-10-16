// Author: Qixiang Yin (yinqx@ihep.ac.cn)
// Purpose: Transfer Neu_Spectrum (E_deposited) to Neu_Spectrum (E_visible).

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TH1D.h"
#include "TH1F.h"

struct Data_Set
{
    double E_depo; // 沉积能量
    double Ratio;  // 可见能量/沉积能量
    double E_vis;  // 可见能量
    double Neu;    // 与沉积能量对应的中微子数
};

int main()
{
    TFile* Depo_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/FY.root");
    TH1F* Depo_Spec = (TH1F*)Depo_File -> Get("h_TotalNeuSpec;1");

    TFile* Ratio_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/JUNOInputs2024_07_22.root");
    TH1D* Ratio_Graph = (TH1D*)Ratio_File -> Get("J22rc0_positronScintNL;1"); //从1.02MeV时开始记录

    double nBin = 10e2;         //沉积能谱的bin数
    double Energy_Bin = 1e-2;

    std::vector<Data_Set> Group;

    for (int i = 1; i <= int(nBin); i++)
    {  
        if (Depo_Spec -> GetBinLowEdge(i) + (1e-7) >= 1.02)
        {
            Group.push_back({0,0,0,0});

            Group.back().E_depo = Depo_Spec -> GetBinCenter(i);
            Group.back().Neu = Depo_Spec -> GetBinContent(i);

            for (int j = 1; j <= Ratio_Graph -> GetNbinsX(); j++)
            {
                // 沉积能谱每个bin的中心能量 == Vis/Depo每个bin的左边界能量
                if (std::abs(Depo_Spec -> GetBinCenter(i) - Ratio_Graph -> GetBinLowEdge(j)) < (1e-6))
                {
                     Group.back().Ratio = Ratio_Graph -> GetBinContent(j);
                }
            }
        }
    }
    
    for (int i = 0; i < Group.size(); i++)
    {
        Group[i].E_vis = Group[i].Ratio * Group[i].E_depo;
    }

    //按E_vis排序
    for (int i = 0; i < Group.size(); i++)
    {
        for(int j = 0; j < Group.size() - i - 1; j++)
        {
            if(Group[j].E_vis > Group[j+1].E_vis)
            {
                std::swap(Group[j], Group[j+1]);
            }
        }
    }

    //对E_vis四舍五入到百分位
    for (int i = 0; i < Group.size(); i++)
    {
        Group[i].E_vis = round(Group[i].E_vis * 100.0) / 100.0;
    }
    
    TH1F* h_NeuSpec_E_vis = new TH1F("h_NeuSpec_E_vis", "", int(11e2 - 92.0), 0.92, 11);
    h_NeuSpec_E_vis -> GetXaxis() -> SetTitle("E_vis[MeV]");

    for (int i = 1; i <= int(11e2 - 92.0); i++)
    {
        h_NeuSpec_E_vis -> SetBinContent(i,0);
    }

    //绘制可见能谱
    for (int i = 0; i < Group.size(); i++)
    {
        for (int j = 1; j <= int(11e2 - 92.0); j++)
        {
            if (std::abs(Group[i].E_vis - h_NeuSpec_E_vis -> GetBinLowEdge(j)) < (1e-6))
            {
                h_NeuSpec_E_vis -> SetBinContent(j,h_NeuSpec_E_vis -> GetBinContent(j) + Group[i].Neu);
            }
        }
    }

    for (int i = 1; i <= int(11e2 - 92.0); i++)
    {
        if (h_NeuSpec_E_vis -> GetBinContent(i) == 0.0 && h_NeuSpec_E_vis -> GetBinContent(i-1) != 0 && h_NeuSpec_E_vis -> GetBinContent(i+1) != 0)
        {
            //若bin没有被填充，则取前后两bin的平均值填充
            h_NeuSpec_E_vis -> SetBinContent(i, 0.5 * (h_NeuSpec_E_vis -> GetBinContent(i-1) + h_NeuSpec_E_vis -> GetBinContent(i+1)));
        }
    }

    TFile* outfile = new TFile("./NeuSpec_E_vis.root", "RECREATE");
    outfile -> cd();

    h_NeuSpec_E_vis -> Write();

    outfile -> Close();


    return 0;
}