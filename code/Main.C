#include <iostream>
#include <string.h>

#include "CaenooReader.h"
#include "CaenooEvent.h"

#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"

void Main()
{
    
	int n;
	//std::cin>>n;
	int i=0;
	int Len;
	int Trig,delt,AverSignalCount=0,delt0,delt1,delt2,delt3,AverSignalCount_All=0;
	double Temp,Temp1,Weight,AverSignalScale=0,AverSignalScale_All=0;
	long double Integral,IntegralMid;
	double Baseline;
	int SyncMax=0;
	int SyncMaxAver=0;
	int SyncStart=0;
	int SyncStartAver=0;
	double SyncCOM=0;
	int SyncCOMAver=0;
	double SyncSlope=0;
	int SyncSlopeAver=0;
	double SyncMax_DifFit=0;
	double SyncMax_SqrFit=0;
	int x_fit_min=-1;
	int x_fit_max=-1;
	int Bad_peak_count=0;

	const std::string FilePath="../../1300/lso.dat";
	// const std::string FilePath="../../1100V/waveform_1100_long.dat";
	// const std::string FilePath="../../1100V/waveform_1100_short.dat";
	double Border_Fract=0.95;
	// double Border_Int_l=16000;
	// double Border_Int_u=18000;
	double Border_Int_l=38000;
	double Border_Int_u=46000;
	int N1=15;//max number of points before Max taken for fit
	int N2=15;//max number of points after Max taken for fit
	double MaxFract=0.2;

	double TopFract=0.1;//Max height dif from max for points taken for fit(in fraction from peak height)
	int Mod=0;//0- no drawing; 1- Draw signal forms of events rejected by minimum Peak/Full integral fraction one by one(press enter to veiw next);2-view the fit of the max +-N points one by one(press enter to veiw next);3-view every event signal form;
	int Mod2=1;//Sync point: 0-peak start(Default); 1-Center of mass; 2-Maximum;3- Half height point of front slope
	int Mod3=2;//Weight - 0-the same normalizing coefficient(Default); 1- 1/(full signal integral); 2- 1/(signal peak integral)
	int Mod4=0;//0-no cut by peak shape;1-cut to remove signals with flat or forked top;
	caenoo::Event e;
    caenoo::Reader r;

	r.SetPathToFile(FilePath);
	AverSignalScale=0;
	AverSignalCount=0;
	r.ReadEvent( e );
	auto pCh = e.GetChannel( 0 );
	Len=pCh->GetLength();
	Trig=e.GetPreTrigger();
	TGraph *Gr = new TGraph(Len);
	r.Reset();

	TH1I *histI=new TH1I("hist1","Full signal integral",250,0,100000);
	TH1I *histIMid=new TH1I("hist2","Signal peak integral",250,0,100000);
	TH2I *histI_Imid=new TH2I("hist3","Full from signal peak integrals;Full integral;Signal peak integral",250,0,100000,250,0,100000);
	TH1F *histP_F=new TH1F("hist4","Peak to full signal fraction",100,0,1);
	TH1I *histSM=new TH1I("hist5","SyncMax",Len,0,Len*2);
	TH1I *histSS=new TH1I("hist6","SyncStart",Len,0,Len*2);
	TH1F *histSC=new TH1F("hist7","SyncCOM",Len,0,Len*2);
	TH1F *histSSl=new TH1F("hist12","SyncSlope",Len,0,Len*2);
	TH1F *histB=new TH1F("hist8","Baseline",100,15630,15650);
	TH1F *histDeltB=new TH1F("hist9","Delta from baseline",200,-50,50);
	TH1F *histSignalAver=new TH1F("histSignalAver","",Len-100,100,Len*2-100);

	TH1F *histSignalAver0=new TH1F("hist10_0","Sync point: peak start",Len-100,100,Len*2-100);
	TH1F *histSignalAver1=new TH1F("hist10_1","Sync point: center of mass",Len-100,100,Len*2-100);
	TH1F *histSignalAver2=new TH1F("hist10_2","Sync point: maximum",Len-100,100,Len*2-100);
	TH1F *histSignalAver3=new TH1F("hist10_3","Sync point: front slope middlepoint",Len-100,100,Len*2-100);
	

	TH1I *histDif=new TH1I("hist11","Differencial signal",200,-100,100);
	TGraph *GrDiff=new TGraph(1+N1+N2);
	GrDiff->SetNameTitle("1111","First derivative");
	TGraph *GrInt=new TGraph(1+N1+N2);
	GrInt->SetNameTitle("111","Integral(in fractions from max)");
	GrInt->SetMinimum(0);
	GrInt->SetMaximum(0.6);
	TGraphErrors *GrAlt=new TGraphErrors(1+N1+N2);
	GrAlt->SetNameTitle("11111","Signal");
	TString HistName="#splitline{Averaged signal ";
	switch (Mod2)
	{
	case 0:
		HistName+="Syncronized by signal start, ";
		break;
	case 1:
		HistName+="Syncronized by signal center of mass, ";
		break;
	case 2:
		HistName+="Syncronized by signal maximum peak, ";
		break;
	case 3:
		HistName+="Syncronized by front slope midpoint, ";
		break;
	default:
		HistName+="Syncronized by signal start, ";
		break;
	}
	switch (Mod3)
	{
	case 0:
		HistName+="equally weighted";
		break;
	case 1:
		HistName+="normalized by full integral";
		break;
	case 2:
		HistName+="normalized by peak integral";
		break;
	default:
		HistName+="equally weighted";
		break;
	}
	switch (Mod4)
	{
	case 0:
		HistName+=" no cut by form";
		break;
	case 1:
		HistName+=" added cut by form";
		break;
	default:
		HistName+=" no cut by form";
		break;
	}
	HistName+=Form("}{signal integral range=(%.0f,%.0f)",Border_Int_l,Border_Int_u);
	

	
	TCanvas* c0=new TCanvas();
	TCanvas* c5=new TCanvas();
	TCanvas* c6=new TCanvas();
	
	while(r.ReadEvent( e ))
    {
		SyncCOM=0;
		Integral=0;
		IntegralMid=0;
		SyncSlope=-1;
		pCh = e.GetChannel( 0 );
		if( pCh == nullptr )
		{
			continue;
		}
		Baseline=0;
		Temp=-1;
		for(int o=0;o<Trig/2;o++)
		{
			Baseline+=2.0*pCh->At(o)/(double)(Trig);
		}
		histB->Fill(Baseline);
		for(int o=0;o<Trig/2;o++)
		{
			histDeltB->Fill(pCh->At(o)-Baseline);
		}
		for(int o=0;o<Len;o++)
		{
			if(pCh->At(o)>Baseline)
			{
				Temp1=o*2;
			}
			Gr->SetPoint(o,o*2,pCh->At(o));
			Integral+=Baseline-pCh->At(o);
			if((pCh->At(o)<Temp)||(Temp==-1))
			{
				SyncMax=o*2;
				Temp=pCh->At(o);
				SyncStart=Temp1;
			}
		}
		for(int o=0;o<Len;o++)
		{
			SyncCOM+=o*2*(Baseline-pCh->At(o))/Integral;
			if(((Baseline-pCh->At(o))/(Baseline-Temp)>0.5)&&(SyncSlope<0))
			{
				SyncSlope=2.0*((o-1)*((Baseline-pCh->At(o))-(Baseline-Temp)*0.5)+(o)*((Baseline-Temp)*0.5-(Baseline-pCh->At(o-1))))/(pCh->At(o-1)-pCh->At(o));
			}
		}
		for(int o=SyncStart/2;(o<SyncStart/2+150)&&(o<Len);o++)
		{
			IntegralMid+=Baseline-pCh->At(o);
		}
		if((double)(IntegralMid)/(double)(Integral)>Border_Fract)
		{
			histSS->Fill(SyncStart);
			histSC->Fill(SyncCOM);
			histSM->Fill(SyncMax);
			histSSl->Fill(SyncSlope);
		}

	}
	SyncCOMAver=histSC->GetMean();
	SyncStartAver=histSS->GetMean();
	SyncMaxAver=histSM->GetMean();
	SyncSlopeAver=histSSl->GetMean();

	Int_t FR_1,FR_2,FR_3;

	r.Reset();
    TLine a(0,0,0,0);
	TLine b(0,0,0,0);
	TLine c(0,0,0,0);
	TLine d(0,0,0,0);
	
	double NoiseErr=histDeltB->GetStdDev();
	a.SetLineColor(kBlue);
	while(r.ReadEvent( e ))
    {
		SyncCOM=0;
		Integral=0;
		IntegralMid=0;
		SyncSlope=-1;
		pCh = e.GetChannel( 0 );
		if( pCh == nullptr )
		{
			continue;
		}
		
		Baseline=0;
		Temp=-1;
		
		for(int o=0;o<Trig/2;o++)
		{
			Baseline+=2.0*pCh->At(o)/(double)(Trig);
		}

		for(int o=0;o<Len;o++)
		{
			if(pCh->At(o)>Baseline)
			{
				Temp1=o*2;
			}
			Gr->SetPoint(o,o*2,pCh->At(o));
			Integral+=Baseline-pCh->At(o);
			if((pCh->At(o)<Temp)||(Temp==-1))
			{
				SyncMax=o*2;
				Temp=pCh->At(o);
				SyncStart=Temp1;
			}
		}
		for(int o=SyncStart/2;(o<SyncStart/2+150)&&(o<Len);o++)
		{
			IntegralMid+=Baseline-pCh->At(o);
		}		
		for(int o=0;o<Len;o++)
		{
			SyncCOM+=o*2*(Baseline-pCh->At(o))/Integral;
			if(((Baseline-pCh->At(o))/(Baseline-Temp)>0.5)&&(SyncSlope<0))
			{
				SyncSlope=2.0*((o-1)*((Baseline-pCh->At(o))-(Baseline-Temp)*0.5)+(o)*((Baseline-Temp)*0.5-(Baseline-pCh->At(o-1))))/(pCh->At(o-1)-pCh->At(o));
			}
		}
		
		histI->Fill(Integral);
		histIMid->Fill(IntegralMid);
		histI_Imid->Fill(Integral,IntegralMid);
		histP_F->Fill((double)(IntegralMid)/(double)(Integral));
		if(((double)(IntegralMid)/(double)(Integral)<Border_Fract)&&(Mod==1))
		{
			c0->Clear();
			c0->cd();
			Gr->Draw();
			c0->Update();
			std::cout<<"event #"<<i<<"\nBaseline= "<<Baseline<<"\nSyncMax= "<<SyncMax<<"\nSyncCOM= "<<SyncCOM<<"\nSyncStart= "<<SyncStart<<"\nSyncSlope= "<<SyncSlope<<"\nPeak Fraction = "<<(double)(IntegralMid)/(double)(Integral)<<"\n";
			std::cin.ignore();
		}
		if(Mod==3)
		{
			c0->Clear();
			c0->cd();
			Gr->Draw();
			c0->Update();
			std::cout<<"event #"<<i<<"\nBaseline= "<<Baseline<<"\nSyncMax= "<<SyncMax<<"\nSyncCOM= "<<SyncCOM<<"\nSyncStart= "<<SyncStart<<"\nSyncSlope= "<<SyncSlope<<"\nPeak Fraction = "<<(double)(IntegralMid)/(double)(Integral)<<"\n";
			std::cin.ignore();
		}
		
		x_fit_min=-1;
		for(int o=(SyncMax/2)-N1;o<=(SyncMax/2)+N2;o++)
		{
			GrDiff->SetPoint(o-((SyncMax/2)-N1),o*2,pCh->At(o+1)-pCh->At(o-1));
			GrAlt->SetPoint(o-((SyncMax/2)-N1),o*2,pCh->At(o));
			GrAlt->SetPointError(o-((SyncMax/2)-N1),0,sqrt(fabs(Baseline-pCh->At(o))));
			if(pCh->At(o)<Temp+(Baseline-Temp)*TopFract)
			{
				if(x_fit_min<0)
				{
					x_fit_min=o*2;
				}
				x_fit_max=o*2;
			}
		}
		Temp1=0;
		for(int o=0;o<Len;o++)
		{
			Temp1+=(Baseline-pCh->At(o))/Integral;
			if((o>=(SyncMax/2)-N1)&&(o<=(SyncMax/2)+N2))
			{
				GrInt->SetPoint(o-((SyncMax/2)-N1),o,Temp1);
			}
		}
		if(x_fit_max-x_fit_min<4)
		{
			x_fit_max+=2;
			x_fit_min-=2;
		}
		if(Mod==2)
		{
			c0->Clear();
			c0->cd();
			FR_1=GrAlt->Fit("pol2","Q","",x_fit_min,x_fit_max);
			GrAlt->Draw("LAP");
			c0->Update();
			a.DrawLine(SyncMax-2*N1,Temp+(Baseline-Temp)*TopFract,SyncMax+2*N2,Temp+(Baseline-Temp)*TopFract);
			c0->Update();
			c5->Clear();
			c5->cd();
			FR_2=GrDiff->Fit("pol1","Q","",x_fit_min,x_fit_max);
			GrDiff->Draw();
			c5->Update();

			c6->Clear();
			c6->cd();
			FR_3=GrInt->Fit("pol1","Q","",x_fit_min/2,x_fit_max/2);
			GrInt->Draw();
			b.DrawLine(SyncMax/2,0,SyncMax/2,0.6);
			c6->Update();
			c0->cd();
			c.DrawLine(SyncMax-2*N1,Baseline-GrInt->GetFunction("pol1")->GetParameter(1)*Integral,SyncMax+2*N2,Baseline-GrInt->GetFunction("pol1")->GetParameter(1)*Integral);
			d.DrawLine(2*(MaxFract-GrInt->GetFunction("pol1")->GetParameter(0))/GrInt->GetFunction("pol1")->GetParameter(1),0,2*(MaxFract-GrInt->GetFunction("pol1")->GetParameter(0))/GrInt->GetFunction("pol1")->GetParameter(1),17000);
			c0->Update();
		}
		else
		{
			FR_1=GrAlt->Fit("pol2","Q","",x_fit_min,x_fit_max);
			FR_2=GrDiff->Fit("pol1","Q","",x_fit_min,x_fit_max);
			FR_3=GrInt->Fit("pol1","Q","",x_fit_min/2,x_fit_max/2);
		}
		SyncMax_DifFit=-1;
		SyncMax_SqrFit=-1;
		if(FR_1==0)
		{
			SyncMax_SqrFit=-0.5*(GrAlt->GetFunction("pol2")->GetParameter(1)/GrAlt->GetFunction("pol2")->GetParameter(2));
		}
		if(FR_2==0)
		{
			SyncMax_DifFit=-1.0*(GrDiff->GetFunction("pol1")->GetParameter(0)/GrDiff->GetFunction("pol1")->GetParameter(1));
		}
		
		if(Mod==2)
		{
			std::cout<<"max="<<SyncMax<<" differential linear fit max="<<SyncMax_DifFit<<" parbola fit max "<<SyncMax_SqrFit<<" Chi^2/N = "<<(GrAlt->GetFunction("pol2")->GetChisquare()/((x_fit_max-x_fit_min)/2))<<" range "<<x_fit_min<<":"<<x_fit_max<<"\n";
			std::cin.ignore();
		}
		if((Mod4==1)&&((GrAlt->GetFunction("pol2")->GetParameter(2)<0)||(SyncMax_SqrFit<x_fit_min-2)||(SyncMax_SqrFit>x_fit_max+2)||((GrAlt->GetFunction("pol2")->GetChisquare())/((x_fit_max-x_fit_min)/2)>1)))
		{
			Bad_peak_count++;
			continue;
		}
		if(((double)(IntegralMid)/(double)(Integral)>Border_Fract)&&((IntegralMid<Border_Int_u)&&(IntegralMid>Border_Int_l)))
		{
			
			switch (Mod2)
			{
			case 0:
				delt=SyncStartAver-SyncStart;
				break;
			case 1:
				delt=SyncCOMAver-SyncCOM;
				break;
			case 2:
				delt=SyncMaxAver-SyncMax;
				break;
			case 3:
				delt=SyncSlopeAver-SyncSlope;
			default:
				delt=SyncStart-SyncStart;
				break;
			}
			

			delt0=SyncStartAver-SyncStart;
			delt1=SyncCOMAver-SyncCOM;
			delt2=SyncMaxAver-SyncMax;
			delt3=SyncSlopeAver-SyncSlope+3;

			switch (Mod3)
			{
			case 0:
				Weight=Len*2/(Border_Int_l+Border_Int_u);
				break;
			case 1:
				Weight=Len*2/Integral;
				break;
			case 2:
				Weight=Len*2/IntegralMid;
				break;
			default:
				Weight=Len*2/(Border_Int_l+Border_Int_u);
				break;
			}
			if(abs(delt)<100)
			{
				for(int o=50;o<Len-50;o++)
				{
					histSignalAver->Fill(o*2,(Baseline-pCh->At((int)((2*o-delt)/2)))*Weight);
					//histSignalAver->AddBinContent(histSignalAver->FindBin(o*2),(Baseline-pCh->At((int)((2*o-delt)/2)))*Weight);

				}
				AverSignalScale+=Weight;
				AverSignalCount++;
			}

			if((abs(delt0)<100)&&(abs(delt2)<100)&&(abs(delt3)<100))
			{
				for(int o=50;o<Len-50;o++)
				{
					//histSignalAver->Fill(o*2,(Baseline-pCh->At((int)((2*o-delt)/2)))*Weight);
					histSignalAver0->AddBinContent(histSignalAver0->FindBin(o*2),(Baseline-pCh->At((int)((2*o-delt0)/2)))*Weight);
					//histSignalAver1->AddBinContent(histSignalAver1->FindBin(o*2),(Baseline-pCh->At((int)((2*o-delt1)/2)))*Weight);
					histSignalAver2->AddBinContent(histSignalAver2->FindBin(o*2),(Baseline-pCh->At((int)((2*o-delt2)/2)))*Weight);
					histSignalAver3->AddBinContent(histSignalAver3->FindBin(o*2),(Baseline-pCh->At((int)((2*o-delt3)/2)))*Weight);
				}
				AverSignalScale_All+=Weight;
				AverSignalCount_All++;
			}

		}
		i++;
    }

	TCanvas *c1= new TCanvas();
	c1->Divide(2,2);
	c1->cd(1);
	histIMid->Draw();
	gPad->Update();
	Double_t line_max_y=0,dud;
	gPad->GetRangeAxis(dud,dud,dud,line_max_y);
	a.SetLineColor(kRed);
    b.SetLineColor(kRed);
	a.DrawLine(Border_Int_l,0,Border_Int_l,line_max_y);
	b.DrawLine(Border_Int_u,0,Border_Int_u,line_max_y);
	gPad->Update();
	c1->cd(2);
	histI_Imid->Draw();
	c1->cd(3);
	histP_F->Draw();
	c1->cd(4);
	histI->Draw();
	gPad->Update();
	
	TCanvas *c2= new TCanvas();
	c2->Divide(2,2);
	c2->cd(3);
	histSM->Draw();
	c2->cd(1);
	histSS->Draw();
	c2->cd(2);
	histSC->Draw();
	c2->cd(4);
	histSSl->Draw();
	TCanvas *c3=new TCanvas();
	c3->Divide(2,1);
	c3->cd(2);
	histB->Draw();
	c3->cd(1);
	TF1 Ftemp("Ftemp","gaus(0)+gaus(3)",-50,50);
	Ftemp.SetParameters(80000,-5,6,80000,5,6);
	histDeltB->Fit("Ftemp","M","",-50,50);
	histDeltB->Draw();
	HistName+=Form(", number of events averaged=%d};nanoseconds",AverSignalCount);
	histSignalAver->SetTitle(HistName);
	TCanvas *c4=new TCanvas();
	c4->cd();
	histSignalAver->Scale(1.0/histSignalAver->GetMaximum());
	histSignalAver->Draw();
	std::cout<<"\n Events rejected for bad form: "<<Bad_peak_count;
	std::cout<<"\nScale= "<<AverSignalScale<<"\nCount= "<<AverSignalCount<<"\n";
	TCanvas *c7=new TCanvas();
	c7->cd();
	histSignalAver2->Draw("SAME");
	histSignalAver0->Draw("SAME");
	//histSignalAver1->Draw("PLC SAME");
	histSignalAver3->Draw("SAME");
	c7->BuildLegend();
	std::cout<<"\nScale= "<<AverSignalScale_All<<"\nCount= "<<AverSignalCount_All<<"\n";
	TString Filename=Form("result%d%d.root",Mod2,Mod4);
	TFile *f1 = new TFile(Filename,"RECREATE");
	histSignalAver->SetDirectory(f1);
	//histSignalAver->Write("histSignalAver");
	f1->Write();
	f1->Close();
}