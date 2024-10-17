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
    
    TH1F* h_Gaussian_Dis = new TH1F("h_Gaussian_Dis", "", 10e2, 0, 10);
    h_Gaussian_Dis -> GetXaxis() -> SetTitle("E_rec[MeV]");

    for (int i = 1; i <= 10e2; i++)
    {
        h_Gaussian_Dis-> SetBinContent(i,0);
    }
    
    // Gaussian(E_rec, E_vis, Sigma_FWHM)
    double E_vis = 5.0;
    double Sigma_FWHM = GetEnergyResolution(Energy_Resolution_Graph, E_vis) * E_vis;

    for (int i = 1; i <= 10e2; i++)
    {
        double E_rec = 0.5 * Energy_Bin + (i-1) * Energy_Bin;
        h_Gaussian_Dis -> SetBinContent(i, h_Gaussian_Dis -> GetBinContent(i) + GaussianFunc(E_rec, E_vis, Sigma_FWHM));
    }

    TFile* outfile = new TFile("./Gaussian_Dis.root", "RECREATE");
    outfile -> cd();

    h_Gaussian_Dis -> Write();

    outfile -> Close();

    return 0;
}
