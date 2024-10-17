// Author: Qixiang Yin (yinqx@ihep.ac.cn)
// Purpose: Gamma(Reconstruct Energy, Visible Energy, FWHM)

#include <iostream>
#include <string>
#include <cmath>
#include "TGraph.h"
#include "TFile.h"
#include "TH1F.h"

double GetEnergyResolution(TGraph* Graph, double Energy)
{
    double Energy_Resolution;
    
    struct Point
    {
        double x;
        double y;
    };

    std::vector<Point> Graph_Points;
    
    int nPoints = Graph -> GetN();

    for (int i = 0; i < nPoints; i++)
    {
        double x, y;
        Graph -> GetPoint(i, x, y);
        Graph_Points.push_back({x, y});
    }

    for (int i = 0; i < nPoints; i++)
    {
        if (Energy > Graph_Points[i].x && Energy < Graph_Points[i+1].x)
        {
            double slope = (Graph_Points[i+1].y - Graph_Points[i].y) / (Graph_Points[i+1].x - Graph_Points[i].x);
            double intercept = Graph_Points[i].y - slope * Graph_Points[i].x;
            Energy_Resolution = slope * Energy + intercept;
        }
        if (Energy == Graph_Points[i].x)
        {
            Energy_Resolution = Graph_Points[i].y;
        }
    }

    return Energy_Resolution;

}

double GaussianFunc(double x, double E, double Sigma)
{
    double norm = 1.0 / (Sigma * std::sqrt(2 * 3.141));
    double exponent = - 0.5 * std::pow((x - E) / Sigma, 2);
    return norm * std::exp(exponent); 
}

int main()
{
    // Input
    TFile* TAO_Input = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/TAOInputs2022_01_06.root");
    TGraph* Energy_Resolution_Graph = (TGraph*)TAO_Input -> Get("TAOEnergyResolutionNoNeutronRecoil;1");

    TFile* Neu_Spec_Vis_Input = TFile::Open("/junofs/users/yinqixiang/work/FineStructure/EnergyResolution/Program/Depo_to_Vis/build/NeuSpec_E_vis.root");
    TH1F* Neu_Spec_Vis_Graph = (TH1F*)Neu_Spec_Vis_Input -> Get("h_NeuSpec_E_vis;1");
    
    double nBin = Neu_Spec_Vis_Graph -> GetNbinsX();        //可见能谱bin数
    double Start = Neu_Spec_Vis_Graph -> GetBinLowEdge(1);  //可见能谱起始能量
    double Energy_Bin = 1e-2;                               //可见能谱bin的宽度

    // Gaussian分布
    std::vector<double> E_vis; //Gaussian分布能量中心值
    std::vector<double> Sigma_FWHM; //Gaussian分布标准差

    for (int i = 1; i <= nBin; i++)
    {
        double Bin_Center = Neu_Spec_Vis_Graph -> GetBinLowEdge(i) + 0.5 * Energy_Bin;

        double Energy_Resolution_x, Energy_Resolution_y;
        Energy_Resolution_Graph -> GetPoint(0, Energy_Resolution_x, Energy_Resolution_y);
        
        if (Bin_Center > Energy_Resolution_x) //如果Bin_Center大于Energy_Resolution_Graph的第一个点的横坐标
        {
            E_vis.push_back(Bin_Center);
        }
    }

    for (int i = 1; i <= E_vis.size(); i++)
    {
        Sigma_FWHM.push_back(GetEnergyResolution(Energy_Resolution_Graph, E_vis[i-1]) * E_vis[i-1]);
    }

    //Gaussian分布直方图
    TH1F* h_Gaussian_Dis[E_vis.size()];

    for (int i = 0; i < E_vis.size(); i++)
    {
        std::string Gaussian_Graph_Name = "h_Gaussian_Dis_" + std::to_string(E_vis[i]);

        h_Gaussian_Dis[i] = new TH1F(Gaussian_Graph_Name.c_str(), "", int(nBin), Start, Start + nBin * Energy_Bin); //Gaussian分布范围与可见能谱保持一致
        h_Gaussian_Dis[i] -> GetXaxis() -> SetTitle("E_rec[MeV]");
        h_Gaussian_Dis[i]-> SetBinContent(i,0);
    }
    
    // Gaussian(E_rec, E_vis, Sigma_FWHM)
    for (int i = 0; i < E_vis.size(); i++)
    {
        for (int j = 0; j < int(nBin); j++)
        {
            //以可见能谱每个bin的中心能量作为Gaussian分布的重建能量E_rec
            double E_rec = Start + 0.5 * Energy_Bin + j * Energy_Bin;
            
            h_Gaussian_Dis[i] -> SetBinContent(j+1, h_Gaussian_Dis[i] -> GetBinContent(j+1) + GaussianFunc(E_rec, E_vis[i], Sigma_FWHM[i]));
        }
    }

    TFile* outfile = new TFile("./Gaussian_Dis.root", "RECREATE");
    outfile -> cd();

    for (int i = 0; i < E_vis.size(); i++)
    {
        h_Gaussian_Dis[i] -> Write();
    }

    outfile -> Close();

    return 0;
}
