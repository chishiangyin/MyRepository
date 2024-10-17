// Author: Qixiang Yin (yinqx@ihep.ac.cn)
// Purpose: 卷积能量分辨率到可见能谱，绘制重建能谱

#include <iostream>
#include <vector>
#include <string>
#include <TKey.h>
#include <TList.h>
#include <TClass.h>
#include "TFile.h"
#include "TH1F.h"
#include "TGraph.h"

int main()
{
    //能量分辨率曲线
    TFile* TAO_Input = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/TAOInputs2022_01_06.root");
    TGraph* Energy_Resolution_Graph = (TGraph*)TAO_Input -> Get("TAOEnergyResolutionNoNeutronRecoil;1");

    //能量分辨率曲线的第一个坐标点
    double Energy_Resolution_x, Energy_Resolution_y;
    Energy_Resolution_Graph -> GetPoint(0, Energy_Resolution_x, Energy_Resolution_y);
    
    //可见能谱
    TFile* Visible_Spec_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/Program/Depo_to_Vis/build/NeuSpec_E_vis.root");
    TH1F* Visible_Spec_TEMP = (TH1F*)Visible_Spec_File -> Get("h_NeuSpec_E_vis;1");

    //E_vis不同的Gaussian分布的集合
    TFile* Gaussian_Resolution_Matrix_File = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/Program/Gaussian_Dis/build/Gaussian_Dis.root");
    std::vector<TH1F*> Gaussian_Resolution_Matrix;

    TList *keyList = Gaussian_Resolution_Matrix_File -> GetListOfKeys();
    TIter next(keyList);
    TKey *key;

    while ((key = (TKey*)next())) 
    {
        if (key -> GetClassName() == std::string("TH1F")) 
        {
            TH1F *hist = (TH1F*)key->ReadObj();
            Gaussian_Resolution_Matrix.push_back(hist);
        }
    }

    double nBin = Visible_Spec_TEMP -> GetNbinsX();   //可见能谱bin数
    double Energy_Bin = 1e-2;                         //可见能谱bin的宽度
    int nBin_Rec = 0;                                 //重建能谱bin数
    double Start;                                     //重建能谱起始能量
    int Flag_1 = 0;

    for (int i = 1; i <= int(nBin); i++)
    {
        if (Visible_Spec_TEMP -> GetBinLowEdge(i) + 0.5 * Energy_Bin > Energy_Resolution_x) //如果bin的中心值大于能量分辨率曲线的第一个点的横坐标
        {
            nBin_Rec++;

            if (Flag_1 == 0)
            {
                Start = Visible_Spec_TEMP -> GetBinLowEdge(i);
                Flag_1 = 1;
            }
        }
    }

    //将原始可见能谱转换为和重建能谱范围一致的能谱
    TH1F* Visible_Spec = new TH1F("h_Neu_Spec_E_vis_Range", "", int(nBin_Rec), Start, Start + double(nBin_Rec) * Energy_Bin);
    
    int Flag_2 = 1;

    for (int i = 1; i <= int(nBin); i++)
    {
        if (Visible_Spec_TEMP -> GetBinLowEdge(i) >= Start)
        {
            Visible_Spec -> SetBinContent(Flag_2, Visible_Spec_TEMP -> GetBinContent(i));
            Flag_2++;
        }
    }

    //绘制重建能谱
    TH1F* Reconstruct_Spec = new TH1F("h_Neu_Spec_E_rec", "", int(nBin_Rec), Start, Start + double(nBin_Rec) * Energy_Bin);
    Reconstruct_Spec -> GetXaxis() -> SetTitle("E_rec[MeV]");

    for (int i = 1; i <= int(nBin_Rec); i++)
    {
        Reconstruct_Spec -> SetBinContent(i, 0);
    }

    for (int i = 0; i < Gaussian_Resolution_Matrix.size(); i++) //循环重建能量
    {
        double Neu_Rec = 0.0;

        for (int j = 0; j < int(nBin_Rec); j++) //循环可见能量
        {
            Neu_Rec += Energy_Bin * Visible_Spec -> GetBinContent(j+1) * Gaussian_Resolution_Matrix[j] -> GetBinContent(i+1); //卷积积分
        }

        Reconstruct_Spec -> SetBinContent(i+1, Neu_Rec);
    }

    TFile* outfile = new TFile("./Convolution_Spec.root", "RECREATE");
    outfile -> cd();

    Reconstruct_Spec -> Write();

    outfile -> Close();

    return 0;
}