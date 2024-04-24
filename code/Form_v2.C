//Version with tau_rise fit
#include <iostream>
#include <string.h>

#include "CaenooReader.h"
#include "CaenooEvent.h"

#include "TGraph.h"
#include "TH1.h"

#include "TCanvas.h"
#include "TStyle.h"

#include "TPad.h"

#include "TF1.h"
#include "TF1Convolution.h"
#include "TFile.h"
#include "TMath.h"

Double_t SignalFuction(Double_t* x, Double_t* p)
{
    Double_t A = p[0];
    Double_t alpha=p[1];
    Double_t tau_r=p[2];
    Double_t tau1=p[3];
    Double_t tau2=p[4];
    Double_t x0=p[5];
    Double_t xx = x[0];
    Double_t Result=0;
    if(xx>x0)
    {
        Result+=A*alpha*TMath::Exp(-(xx-x0)/tau1);
        Result+=A*(1-alpha)*TMath::Exp(-(xx-x0)/tau2);
    }
    return Result*(1-TMath::Exp(-(xx-x0)/tau_r));
}
Double_t Gaus_f(Double_t* x, Double_t* p)
{
    
    Double_t C=p[0];
    Double_t Mean=p[1];
    Double_t Sigma=p[2];
    Double_t xx=x[0];
    Double_t Result=C*TMath::Exp(-TMath::Power((xx-Mean)/Sigma,2)/2)/TMath::Sqrt(2*TMath::Pi())/Sigma;
    return Result;
}
void Form_v2()
{

    /*
    double C=0.012;
    double R=1;
    double I1=10;
    double I2=0;
    double tau1=5;
    double tau2=100;
    double x0=1860;
    TF1 *Signal=new TF1("Signal","((x-[x0])>0?1:0)*(([C]*[I1]/([C]/[R]-1/[tau1]))*(exp(-(x-[x0])/[tau1])-exp(-(x-[x0])*[C]/[R]))+([C]*[I2]/([C]/[R]-1/[tau2]))*(exp(-(x-[x0])/[tau2])-exp(-(x-[x0])*[C]/[R])))",1800,2300);
    Signal->SetParameter("C",C);
    Signal->SetParameter("tau1",tau1);
    Signal->SetParameter("I1",I1);
    Signal->SetParameter("tau2",tau2);
    Signal->SetParameter("I2",I1);
    Signal->SetParameter("R",R);
    Signal->SetParameter("x0",x0);
    */

    int Mod2=2;
    int Mod4=1;

    double A=2;
    double alpha=0.95;
    double tau1=40;
    double tau2=100;
    double tau_r=11;
    double x0=1889;

    double C=1;
    double Sigma=6;

    TF1 *Signal=new TF1("Signal",SignalFuction,-300,4000,6);
    Signal->SetNpx(4301);
    TF1 *Gaus_C=new TF1("Gaus_c",Gaus_f,-50,50,3);
    Signal->SetNpx(101);
    Signal->SetParameter(0,A);
    Signal->SetParameter(1,alpha);
    Signal->SetParameter(2,tau_r);
    Signal->SetParameter(3,tau1);
    Signal->SetParameter(4,tau2);
    Signal->SetParameter(5,x0);
    Gaus_C->SetParameter(0,C);
    Gaus_C->SetParameter(1,0);
    Gaus_C->SetParameter(2,Sigma);
    //Signal->Draw();
    //TH1F *histInitSig=new TH1F("Hist init signal","initial signal",10000,-500,500);
    TF1Convolution *Signal_final=new TF1Convolution(Signal,Gaus_C,-300,4000);
    Signal_final->SetNofPointsFFT(10000);
    //Signal_final->SetRange(-300,4000);

    //std::cout<<Signal_final->GetNpar();
    
    TF1 *Signal_final_2= new TF1("S_F_2",Signal_final,-300,4000,Signal_final->GetNpar());
    Signal_final_2->SetParNames("A","#alpha","#tau_{rise}","#tau_{1}","#tau_{2}","x_{0}","gaus-C","gaus-mean","#sigma_{g}");    
    Signal_final_2->SetParameter(0,A);
    Signal_final_2->SetParameter(1,alpha);
    Signal_final_2->SetParLimits(1,0.0,1.0);
    Signal_final_2->SetParameter(2,tau_r);
    Signal_final_2->SetParameter(3,tau1);
    Signal_final_2->SetParameter(4,tau2);
    //Signal_final_2->SetParLimits(4,500,10000);
    Signal_final_2->SetParameter(5,x0);
    Signal_final_2->SetParameter(8,Sigma);
    Signal_final_2->SetParLimits(8,0.5,100);
//    Signal_final_2->SetParameter(6,C);
//    Signal_final_2->SetParameter(7,0);
    Signal_final_2->FixParameter(6,C);
    Signal_final_2->FixParameter(7,0);
    
    TString Filename=Form("result%d%d.root",Mod2,Mod4);
    TFile *f1= new TFile(Filename,"READ");
    TH1F *histSignal=new TH1F();
    f1->GetObject("histSignalAver",histSignal);
    TCanvas *c1=new TCanvas();

    Gaus_C->Draw();
    TCanvas *c2=new TCanvas();
    Signal->SetNpx(4300);
    Signal->Draw();
    // TCanvas *c3=new TCanvas();
    // Signal_final_2->SetNpx(4300);
    // Signal_final_2->Draw();

    TCanvas *c4=new TCanvas();
    histSignal->Draw();
    Signal_final_2->SetNpx(4300);
    Signal_final_2->Draw("same");
    histSignal->Fit(Signal_final_2,"","",1850,2300);
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);
    TString FilenameOut="IntegralFits/Long_1100";
    switch (Mod2)
	{
	case 0:
		FilenameOut+="_start";
		break;
	case 1:
		FilenameOut+="_com";
		break;
    case 2:
		FilenameOut+="_max";
		break;
	case 3:
		FilenameOut+="_slope";
		break;
            
	default:
		FilenameOut+="_start";
		break;
	}
    switch (Mod4)
	{
	case 0:
		FilenameOut+="_nocut";
		break;
	case 1:
		FilenameOut+="_cut";
		break;
	default:
		FilenameOut+="_nocut";
		break;
	}
    FilenameOut+="_v2.root";
	TFile *f2 = new TFile(FilenameOut,"RECREATE");
	histSignal->SetDirectory(f2);
    f2->Write();
	//f2->Close();
}