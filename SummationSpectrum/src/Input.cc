#include "Input.hh"
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <sstream>
#include <fstream>
Input::Input()
{

}
Input::~Input()
{
    IsoActivityList.clear();
    ENSDFBetaDecayData.clear();
    IsoFYList.clear();
}
std::vector<Input::ISOFYStruct> Input::GetFYList()
{
    return IsoFYList;
}
void Input::SetFYList(Input::ISOFYStruct fy)
{
    IsoFYList.push_back(fy);
}
std::vector<Input::ISOStruct> Input::GetIsoList()
{
    return IsoActivityList;
}
void Input::SetIsoList(Input::ISOStruct iso)
{
    IsoActivityList.push_back(iso);
}
std::vector<Input::ENSDFStruct> Input::GetENSDFData()
{
    return ENSDFBetaDecayData;
}
void Input::SetENSDFData(Input::ENSDFStruct ENSDF)
{
    ENSDFBetaDecayData.push_back(ENSDF);
}
void Input::SetIsoFileName(std::string filename)
{
    isofilename = filename;
}
void Input::SetENSDFFileName(std::string filename)
{
    ENSDFfilename = filename;
}
void Input::SetFYFileName(std::string filename)
{
    FYfilename = filename;
}
void Input::ReadENSDFInputFile()
{
    std::ifstream chanENSDF(ENSDFfilename.c_str());
    std::string lineENSDF;
    Input::ENSDFStruct ensdf;
    std::string ENSDFisoname;
    std::string stringA;
    std::string stringZ;
    std::string stringQvalue;
    std::string stringBranch;
    std::string stringEnergyLevel;
    std::string stringJpi;
    std::string stringParentEnergyLevel;
    std::string stringParentJpi;
    int j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12,j13,j14,j15,j16,j17;
    while(std::getline(chanENSDF, lineENSDF))
    {
        if(lineENSDF.size()>=2)
        {
            if(lineENSDF.at(0)=='B' && lineENSDF.at(1)=='-')
            {
                ENSDFisoname.clear();
                stringA.clear();
                stringZ.clear();
                stringQvalue.clear();
                stringBranch.clear();
                stringEnergyLevel.clear();
                stringJpi.clear();
                stringParentEnergyLevel.clear();
                stringParentJpi.clear();
                for(j1=4;j1<lineENSDF.size();j1++)
                {
                    if(lineENSDF[j1]!=' ')
                    {
                        ENSDFisoname+=lineENSDF[j1];
                    }
                    else
                    {   
                        break;
                    }
                }
                for(j2=j1;j2<lineENSDF.size();j2++)
                {
                    if(lineENSDF[j2]=='A'&&lineENSDF[j2+1]==':')
                    {
                        for(j3=j2+3;j3<lineENSDF.size();j3++)
                        {
                            if(lineENSDF[j3]!=' ')
                            {
                                stringA+=lineENSDF[j3];
                            }
                            else
                            {
                                break;
                            }
                        }
                        break; 
                       
                    }
                }
                for(j4=j3;j4<lineENSDF.size();j4++)
                {
                    if(lineENSDF[j4]=='Z'&&lineENSDF[j4+1]==':')
                    {
                        for(j5=j4+3;j5<lineENSDF.size();j5++)
                        {
                            if(lineENSDF[j5]!=' ')
                            {
                                stringZ+=lineENSDF[j5];
                            }
                            else
                            {
                                break;
                            }
                        }
                        break;
                    }
                }
                for(j6=j5;j6<lineENSDF.size();j6++)
                {
                    if(lineENSDF[j6]=='Q'&&lineENSDF[j6+1]=='-'&&lineENSDF[j6+2]=='v'&&lineENSDF[j6+3]=='a')
                    {
                        for(j7=j6+9;j7<lineENSDF.size();j7++)
                        {
                            if(lineENSDF[j7]!=' ')
                            {
                            
                                stringQvalue+=lineENSDF[j7];
                            }
                            else
                            {
                                break;
                            }
                        }
                        break;
                    }
                }
                for(j8=j7;j8<lineENSDF.size();j8++)
                {
                    if(lineENSDF[j8]=='B'&&lineENSDF[j8+1]=='r'&&lineENSDF[j8+2]=='a'&&lineENSDF[j8+3]=='n')
                    {
                        for(j9=j8+17;j9<lineENSDF.size();j9++)
                        {
                        
                            if(lineENSDF[j9]!=' ')
                            {   
                                stringBranch+=lineENSDF[j9];
                            }
                            else
                            {
                                break;
                            }
                        }
                        break;
                    }
                }
                for(j10=j9;j10<lineENSDF.size();j10++)
                {
                    if(lineENSDF[j10]=='G'&&lineENSDF[j10+1]=='r'&&lineENSDF[j10+2]=='o'&&lineENSDF[j10+3]=='u')
                    {
                        for(j11=j10+22;j11<lineENSDF.size();j11++)
                        {
                        
                            if(lineENSDF[j11]!=' ')
                            {
                                stringEnergyLevel+=lineENSDF[j11];
                            }
                            else
                            {
                                break;
                            }
                        }
                        break;
                    }
                }
                for(j12=j11;j12<lineENSDF.size();j12++)
                {
                    if(lineENSDF[j12]=='J'&&lineENSDF[j12+1]=='p'&&lineENSDF[j12+2]=='i'&&lineENSDF[j12+3]==':')
                    {
                        for(j13=j12+4;j13<lineENSDF.size();j13++)
                        {
                            if(!(lineENSDF[j13]=='P'&&lineENSDF[j13+1]=='a'&&lineENSDF[j13+2]=='r'))
                            {
                                stringJpi+=lineENSDF[j13];
                            }
                            else
                            {
                                break;
                            }
                        }
                        break;
                    }
                }
                for(j14=j13-10;j14<lineENSDF.size();j14++)
                {
                
                    if(lineENSDF[j14]=='P'&&lineENSDF[j14+1]=='a'&&lineENSDF[j14+2]=='r'&&lineENSDF[j14+3]=='e')
                    {
                    
                        for(j15=j14+10;j15<lineENSDF.size();j15++)
                        {
                        
                            if(isdigit(lineENSDF[j15])||lineENSDF[j15]=='.')
                            {
                                for(j16=j15;j16<lineENSDF.size();j16++)
                                {
                                    if(isdigit(lineENSDF[j16])||lineENSDF[j16]=='.'||lineENSDF[j16]=='e'||lineENSDF[j16]=='E'||lineENSDF[j16]=='-'||lineENSDF[j16]=='+')
                                    //if(ENSDFBetaDecay.at(ii).at(ii16)!=' ')
                                    {
                                        stringParentEnergyLevel+=lineENSDF[j16];
                                    }
                                    else if(lineENSDF[j16]==' ')
                                    {
                                        break;
                                    }
                                }
                                break;
                            }

                            if(lineENSDF[j15]=='X'||lineENSDF[j15]=='x')
                            {
                                stringParentEnergyLevel+="0";
                                break;
                            }
                        }
                        break;
                    }
                }
                for(j17=j16;j17<lineENSDF.size();j17++)
                {
                    stringParentJpi+=lineENSDF[j17];
                }
                ensdf.isoname = ENSDFisoname;
                ensdf.A = stoi(stringA);
                ensdf.Z = stoi(stringZ);
                ensdf.Qvalue = stod(stringQvalue)/1000.0;
                ensdf.BranchStrength = stod(stringBranch)*0.01;
                ensdf.EnergyLevel = stod(stringEnergyLevel)/1000.0;
                ensdf.ParentEnergyLevel = stod(stringParentEnergyLevel)/1000.0;
                int k1,k2,k3;
                for(k1=0;k1<stringJpi.size();k1++)
                {
                    if(stringJpi[k1]==' ')
                    {
                        if(k1>=15 || k1==stringJpi.size()-1)
                        {
                            ensdf.JpiFlag = 1;
                            ensdf.J = "Unknow";
                            ensdf.Pi = 4;
                            break;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else
                    {
                        std::string findstring1 = "TO";
                        std::string findstring2 = "To";
                        std::string findstring3 = "to";
                        std::string findstring4 = ",";
                        std::string stringmemory;
                        stringmemory.clear();
                        int index1 =  stringJpi.find(findstring1,0);
                        int index2 =  stringJpi.find(findstring2,0);
                        int index3 =  stringJpi.find(findstring3,0);
                        int index4 =  stringJpi.find(findstring4,0);
                        if((index1<stringJpi.length())||(stringJpi.length())||(index3<stringJpi.length())||(index4<stringJpi.length()))
                        {
                            ensdf.JpiFlag = 2;
                            ensdf.J = stringJpi;
                            ensdf.Pi = 3;
                            break;
                        }
                        else
                        {
                            for(k2=k1;k2<stringJpi.size();k2++)
                            {
                                if(isdigit(stringJpi[k2]))
                                {
                                    for(k3=k2;k3<stringJpi.size();k3++)
                                    {
                                        if(isdigit(stringJpi[k3])||stringJpi[k3]=='/')
                                        {
                                            stringmemory+=stringJpi[k3];
                                        }
                                        else
                                        {
                                            ensdf.J = stringmemory;
                                            break;
                                        }
                                    }

                                    std::string findstring5 = "+";
                                    std::string findstring6 = "-";
                                    int index5 = stringJpi.find(findstring5,k2);
                                    int index6 = stringJpi.find(findstring6,k2);
                                    if((index5<stringJpi.length())||(index6<stringJpi.length()))
                                    {
                                        ensdf.JpiFlag = 0;
                                        if(index5<stringJpi.length())
                                        {
                                            ensdf.Pi = 1;
                                        }

                                        if(index6<stringJpi.length())
                                        {
                                            ensdf.Pi = -1;
                                        }
                                    }
                                    else
                                    {
                                        ensdf.JpiFlag = 3;
                                        ensdf.Pi = 2;
                                    }
                                    break;
                                }
                            }

                            break;
                        }
                    }
                }

                for(k1=0;k1<stringParentJpi.size();k1++)
                {
                    if(stringParentJpi[k1]==' ')
                    {
                        if(k1>=15 || k1==stringParentJpi.size()-1)
                        {
                            ensdf.ParentJpiFlag = 1;
                            ensdf.ParentJ = "Unknow";
                            ensdf.ParentPi = 4;
                            break;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    else
                    {
                        std::string findstring1 = "TO";
                        std::string findstring2 = "To";
                        std::string findstring3 = "to";
                        std::string findstring4 = ",";
                        std::string stringmemory;
                        stringmemory.clear();
                        int index1 =  stringParentJpi.find(findstring1,0);
                        int index2 =  stringParentJpi.find(findstring2,0);
                        int index3 =  stringParentJpi.find(findstring3,0);
                        int index4 =  stringParentJpi.find(findstring4,0);
                        if((index1<stringParentJpi.length())||(stringParentJpi.length())||(index3<stringParentJpi.length())||(index4<stringParentJpi.length()))
                        {
                            ensdf.ParentJpiFlag = 2;
                            ensdf.ParentJ = stringParentJpi;
                            ensdf.ParentPi = 3;
                            break;
                        }
                        else
                        {
                            for(k2=k1;k2<stringParentJpi.size();k2++)
                            {
                                if(isdigit(stringParentJpi[k2]))
                                {
                                    for(k3=k2;k3<stringParentJpi.size();k3++)
                                    {
                                        if(isdigit(stringParentJpi[k3])||stringParentJpi[k3]=='/')
                                        {
                                            stringmemory+=stringParentJpi[k3];
                                        }
                                        else
                                        {
                                            ensdf.ParentJ = stringmemory;
                                            break;
                                        }
                                    }

                                    std::string findstring5 = "+";
                                    std::string findstring6 = "-";
                                    int index5 = stringParentJpi.find(findstring5,k2);
                                    int index6 = stringParentJpi.find(findstring6,k2);
                                    if((index5<stringParentJpi.length())||(index6<stringParentJpi.length()))
                                    {
                                        ensdf.ParentJpiFlag = 0;
                                        if(index5<stringParentJpi.length())
                                        {
                                            ensdf.ParentPi = 1;
                                        }

                                        if(index6<stringParentJpi.length())
                                        {
                                            ensdf.ParentPi = -1;
                                        }
                                    }
                                    else
                                    {
                                        ensdf.ParentJpiFlag = 3;
                                        ensdf.ParentPi = 2;
                                    }
                                    break;
                                }
                            }

                            break;
                        }
                    }
                }
                SetENSDFData(ensdf);


            }
        }
    }
    
    
}
void Input::ReadIsoInputFile()
{
    std::ifstream chanISO(isofilename.c_str());
    std::string lineISO;
    Input::ISOStruct iso;
    double activity;
    std::string sactivity;
    std::string isoname;
    int i1,i2,i3;
    while(std::getline(chanISO, lineISO))
    {
        if(lineISO.size()>=4)
        {
            if(lineISO[0]=='N'&& lineISO[1]=='a' && lineISO[2]=='m'&& lineISO[3]=='e')
            {
                isoname.clear();
                sactivity.clear();
                for(i1=6;i1<lineISO.size();i1++)
                {
                    if(lineISO[i1]!=' ')
                    {
                        isoname+=lineISO[i1];
                    }
                    else
                    {
                        break;
                    }
                    
                }
                for(i2=i1;i2<lineISO.size();i2++)
                {
                    if(lineISO[i2]=='A'&&lineISO[i2+1]=='c'&&lineISO[i2+2]=='t')
                    {
                        for(i3=i2+10;i3<lineISO.size();i3++)
                        {
                            sactivity+=lineISO[i3];
                        }
                        break;
                    }
                }
                activity = std::stod(sactivity);
                iso.isoname = isoname;
                iso.activity = activity;
                SetIsoList(iso);
            }
        }
        
    }
}

void Input::ReadFYInputFile()
{
    std::ifstream chanFY(FYfilename.c_str());
    int count = 0;
    std::string line;
    Input::ISOFYStruct fy;
    while (std::getline(chanFY, line)) {
        std::vector<std::string> row;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            row.push_back(cell);
        }
        if(count>0)
        {
            
            std::string ParentName;
            std::string DaughterName;
            ParentName = row[4];
            ParentName += row[5];
            DaughterName = row[1];
            DaughterName += row[2];
            for (char &c : ParentName) {
                if (std::islower(c)) {
                    c = std::toupper(c);
                }
            }
            for (char &c : DaughterName) {
                if (std::islower(c)) {
                    c = std::toupper(c);
                }
            }
            if(row[7].size()==0)
            {
                row[7]="0";
            }
            if(row[12].size()==0)
            {
                row[12]="0";
            }
            if(row[8].size()==0)
            {
                row[8]="0";
            }
            if(row[9].size()==0)
            {
                row[9]="0";
            }
            if(row[10].size()==0)
            {
                row[10]="0";
            }
            if(row[11].size()==0)
            {
                row[11]="0";
            }
            fy.Parentisoname = ParentName;
            fy.isoname = DaughterName;
            fy.ParentA = std::stoi(row[4]);
            fy.ParentZ = std::stoi(row[3]);
            fy.A = std::stoi(row[1]);
            fy.Z = std::stoi(row[0]);
            fy.NEnergyLevel = std::stoi(row[6]);
            fy.thermal_fy = std::stod(row[7]);
            fy.thermal_fy_unc= std::stod(row[8]);
            fy.fast_fy= std::stod(row[9]);
            fy.fast_fy_unc= std::stod(row[10]);
            fy.fourteen_mev_fy= std::stod(row[11]);
            fy.fourteen_mev_fy_unc= std::stod(row[12]);
            SetFYList(fy);
        }
        count++;
        row.clear();
    }
}
