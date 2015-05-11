#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>

#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TGaxis.h>
#include <TF1.h>
#include <TStyle.h>

#include "/cmshome/fpreiato/DiJet/test/CMSSW_7_2_1_DiJet/src/CMSDIJET/DijetRootTreeAnalyzer/Compare.h"

using namespace std;

///////////////////////////


int main(int argc, char* argv[]) {
  
  ///////////////// WHAT IT HAVE TO DO
  bool Draw_SmearingFunctions_Bin100_SixCategory = false;
  bool Draw_SmearingFunctions_Bin200_SixCategory = false;
  bool Draw_SmearingFunctions_Bin200_FourCategory = false ;
  bool Draw_SmearingFunctions_Bin200_FourCategory_Alone = false ;
  bool Draw_SmearingFunctions_Bin200_TwoCategory_Example = false;
  ///////////////////////////////////////////////////////////////////////////
    
  if(Draw_SmearingFunctions_Bin100_SixCategory == true){

    TFile file1("SmearingFunctions_RecoGen_GravitonSamples/BinPt100GeV/SixCategory/SmearingFunctions_ResonanceSamplesToQuarkQuark_Bin100.root ");
    TH1D *H_step_pt = (TH1D*)file1.Get("H_step_pt");
    int step_pt = H_step_pt->GetMean();
    TFile file2("SmearingFunctions_RecoGen_QCDSamples/BinPt100GeV/SixCategory/SmearingFunctions_QCDSamples_Bin100.root ");    
    
    int n_bin = 4500 / step_pt ;   
    char HistoName[200];	 
    TH1D *h[100][100];
    TH1D *j[100][100];
    
    TLegend *leg = new TLegend( 0.9, 0.9, 0.6, 0.6);
    leg -> SetHeader("Smearing Functions: ");
    leg -> AddEntry(h[0][0], "Resonance Samples", "L");
    leg -> AddEntry(j[0][0], "QCD Samples", "L");
    
    for( int kk=0; kk<6; kk++){//tipi di partone   
      for(int ii=0; ii<5; ii++){//bin in eta
	for(int jj=0; jj<n_bin ; jj++){//bin in pT	 	 
	  int pt_bin_max = (jj*step_pt)+step_pt;
	  // da cambiare e metterci i nomi dei quark
	  if( kk == 0 )	sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 1 )	sprintf(HistoName,"Histo_RQuarkDOWN_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 2 )	sprintf(HistoName,"Histo_RQuarkUP_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 3 )	sprintf(HistoName,"Histo_RQuarkSTRANGE_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 4 )	sprintf(HistoName,"Histo_RQuarkCHARM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 5 )	sprintf(HistoName,"Histo_RQuarkBOTTOM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  j[ii][jj] = (TH1D*)file2.Get(HistoName);   
	  cout<<"      Graviton->Integral(): "<<       h[ii][jj]->Integral() <<endl;
	  cout<<"      QCD->Integral(): "<<       j[ii][jj]->Integral() <<endl;
	  Normalizer(h[ii][jj]);
	  Normalizer( j[ii][jj]);
	  char output_Dir[200] = "Histogram_SmearingFunctions/BinPt100GeV/SixCategory/";
	  char extension[10] = ".png";
	  strcat( HistoName, extension);
	  cout<<HistoName<<endl;
	  if( j[ii][jj]->Integral() >= h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, j[ii][jj], h[ii][jj], 0., 2.,"R","Normalized", kBlue, kRed);
	  if( j[ii][jj]->Integral() < h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, h[ii][jj], j[ii][jj], 0., 2.,"R","Normalized", kRed, kBlue );

	}//end bin pT
      }//end bin eta
    }//end parton type

    //Response in function of flavour for QCD sample
    TH1D *Response_UP            = (TH1D*)file2.Get("Histo_RQuarkUP_WideJet_GenWideJet_0_1600");   
    TH1D *Response_DOWN      = (TH1D*)file2.Get("Histo_RQuarkDOWN_WideJet_GenWideJet_0_1600");   
    TH1D *Response_STRANGE = (TH1D*)file2.Get("Histo_RQuarkSTRANGE_WideJet_GenWideJet_0_1600");   
    TH1D *Response_CHARM    = (TH1D*)file2.Get("Histo_RQuarkCHARM_WideJet_GenWideJet_0_1600");   
    TH1D *Response_BOTTOM  = (TH1D*)file2.Get("Histo_RQuarkBOTTOM_WideJet_GenWideJet_0_1600");   
    TH1D *Response_GLUON    = (TH1D*)file2.Get("Histo_RGluons_WideJet_GenWideJet_0_1600");   
    
    TCanvas *canvas = new TCanvas("canvas","",800,800);
    Response_GLUON    -> SetTitle(" ");
    Response_GLUON    -> SetStats(0);
    Response_GLUON     ->GetXaxis()->SetRangeUser(0.4, 1.3);
    Response_GLUON    -> SetXTitle("Response");
    Response_GLUON    -> SetYTitle("Normalized to unit");
    Response_GLUON    -> GetYaxis()->SetTitleOffset(1.55);
    Response_UP                -> SetLineColor(kRed);
    Response_DOWN         -> SetLineColor(kBlue);
    Response_STRANGE    -> SetLineColor(kGreen+1);
    Response_CHARM       -> SetLineColor(kViolet);
    Response_BOTTOM     -> SetLineColor(kYellow+1);
    Response_GLUON       -> SetLineColor(kBlack);    
    Response_UP                -> SetMarkerColor(kRed);
    Response_DOWN         -> SetMarkerColor(kBlue);
    Response_STRANGE    -> SetMarkerColor(kGreen+1);
    Response_CHARM       -> SetMarkerColor(kViolet);
    Response_BOTTOM     -> SetMarkerColor(kYellow+1);
    Response_GLUON       -> SetMarkerColor(kBlack);    
    Response_UP               -> SetMarkerStyle(20);
    Response_DOWN        -> SetMarkerStyle(21);
    Response_STRANGE   -> SetMarkerStyle(22);
    Response_CHARM      -> SetMarkerStyle(23);
    Response_BOTTOM    -> SetMarkerStyle(31);
    Response_GLUON       -> SetMarkerStyle(33);   
    TLegend *leg2 = new TLegend( 0.1, 0.9, 0.4, 0.6);
    leg2 -> SetHeader("Parton Flavour:");
    leg2 -> AddEntry(Response_UP, "quark Up", "L");
    leg2 -> AddEntry(Response_DOWN, "quark Down", "L");
    leg2 -> AddEntry(Response_STRANGE, "quark Strange", "L");
    leg2 -> AddEntry(Response_CHARM, "quark Charm", "L");
    leg2 -> AddEntry(Response_BOTTOM, "quark Bottom", "L");
    leg2 -> AddEntry(Response_GLUON, "Gluons", "L");
    Response_GLUON      -> DrawNormalized();
    Response_UP              -> DrawNormalized("same");
    Response_DOWN        -> DrawNormalized("same");
    Response_STRANGE   -> DrawNormalized("same");
    Response_CHARM      -> DrawNormalized("same");
    Response_BOTTOM    -> DrawNormalized("same");
    leg2 -> Draw();
    canvas ->SaveAs("Histogram_SmearingFunctions/BinPt100GeV/SixCategory/Compare_Response.png");
    
  }//end if bool

  //////////////////////////////////////////////////////////////////////////////////

  if(Draw_SmearingFunctions_Bin200_SixCategory == true){

    TFile file1("SmearingFunctions_RecoGen_GravitonSamples/BinPt200GeV/SixCategory/SmearingFunctions_ResonanceSamplesToQuarkQuark_Bin200.root ");
    TH1D *H_step_pt = (TH1D*)file1.Get("H_step_pt");
    int step_pt = H_step_pt->GetMean();
    TFile file2("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/SixCategory/SmearingFunctions_QCDSamples_Bin200.root ");    
    
    int n_bin = 4500 / step_pt ;   
    char HistoName[200];	 
    TH1D *h[100][100];
    TH1D *j[100][100];
    
    TLegend *leg = new TLegend( 0.9, 0.9, 0.6, 0.6);
    leg -> SetHeader("Smearing Functions: ");
    leg -> AddEntry(h[0][0], "Resonance Samples", "L");
    leg -> AddEntry(j[0][0], "QCD Samples", "L");
    
    for( int kk=0; kk<6; kk++){//tipi di partone   
      for(int ii=0; ii<5; ii++){//bin in eta
	for(int jj=0; jj<n_bin ; jj++){//bin in pT	 	 
	  int pt_bin_max = (jj*step_pt)+step_pt;
	  // da cambiare e metterci i nomi dei quark
	  if( kk == 0 )	sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 1 )	sprintf(HistoName,"Histo_RQuarkDOWN_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 2 )	sprintf(HistoName,"Histo_RQuarkUP_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 3 )	sprintf(HistoName,"Histo_RQuarkSTRANGE_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 4 )	sprintf(HistoName,"Histo_RQuarkCHARM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 5 )	sprintf(HistoName,"Histo_RQuarkBOTTOM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  j[ii][jj] = (TH1D*)file2.Get(HistoName);   
	  cout<<"      Graviton->Integral(): "<<       h[ii][jj]->Integral() <<endl;
	  cout<<"      QCD->Integral(): "<<       j[ii][jj]->Integral() <<endl;
	  Normalizer(h[ii][jj]);
	  Normalizer( j[ii][jj]);
	  char output_Dir[200] = "Histogram_SmearingFunctions/BinPt200GeV/SixCategory/";
	  char extension[10] = ".png";
	  strcat( HistoName, extension);
	  cout<<HistoName<<endl;
	  if( j[ii][jj]->Integral() >= h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, j[ii][jj], h[ii][jj], 0., 2.,"R","Normalized", kBlue, kRed);
	  if( j[ii][jj]->Integral() < h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, h[ii][jj], j[ii][jj], 0., 2.,"R","Normalized", kRed, kBlue );

	}//end bin pT
      }//end bin eta
    }//end parton type

    //Response in function of flavour for QCD sample
    TH1D *Response_UP            = (TH1D*)file2.Get("Histo_RQuarkUP_WideJet_GenWideJet_0_1600");   
    TH1D *Response_DOWN      = (TH1D*)file2.Get("Histo_RQuarkDOWN_WideJet_GenWideJet_0_1600");   
    TH1D *Response_STRANGE = (TH1D*)file2.Get("Histo_RQuarkSTRANGE_WideJet_GenWideJet_0_1600");   
    TH1D *Response_CHARM    = (TH1D*)file2.Get("Histo_RQuarkCHARM_WideJet_GenWideJet_0_1600");   
    TH1D *Response_BOTTOM  = (TH1D*)file2.Get("Histo_RQuarkBOTTOM_WideJet_GenWideJet_0_1600");   
    TH1D *Response_GLUON    = (TH1D*)file2.Get("Histo_RGluons_WideJet_GenWideJet_0_1600");   
    
    TCanvas *canvas = new TCanvas("canvas","",800,800);
    Response_GLUON    -> SetTitle(" ");
    Response_GLUON    -> SetStats(0);
    Response_GLUON     ->GetXaxis()->SetRangeUser(0.4, 1.3);
    Response_GLUON    -> SetXTitle("Response");
    Response_GLUON    -> SetYTitle("Normalized to unit");
    Response_GLUON    -> GetYaxis()->SetTitleOffset(1.55);
    Response_UP                -> SetLineColor(kRed);
    Response_DOWN         -> SetLineColor(kBlue);
    Response_STRANGE    -> SetLineColor(kGreen+1);
    Response_CHARM       -> SetLineColor(kViolet);
    Response_BOTTOM     -> SetLineColor(kYellow+1);
    Response_GLUON       -> SetLineColor(kBlack);    
    Response_UP                -> SetMarkerColor(kRed);
    Response_DOWN         -> SetMarkerColor(kBlue);
    Response_STRANGE    -> SetMarkerColor(kGreen+1);
    Response_CHARM       -> SetMarkerColor(kViolet);
    Response_BOTTOM     -> SetMarkerColor(kYellow+1);
    Response_GLUON       -> SetMarkerColor(kBlack);    
    Response_UP               -> SetMarkerStyle(20);
    Response_DOWN        -> SetMarkerStyle(21);
    Response_STRANGE   -> SetMarkerStyle(22);
    Response_CHARM      -> SetMarkerStyle(23);
    Response_BOTTOM    -> SetMarkerStyle(31);
    Response_GLUON       -> SetMarkerStyle(33);   
    TLegend *leg2 = new TLegend( 0.1, 0.9, 0.4, 0.6);
    leg2 -> SetHeader("Parton Flavour:");
    leg2 -> AddEntry(Response_UP, "quark Up", "L");
    leg2 -> AddEntry(Response_DOWN, "quark Down", "L");
    leg2 -> AddEntry(Response_STRANGE, "quark Strange", "L");
    leg2 -> AddEntry(Response_CHARM, "quark Charm", "L");
    leg2 -> AddEntry(Response_BOTTOM, "quark Bottom", "L");
    leg2 -> AddEntry(Response_GLUON, "Gluons", "L");
    Response_GLUON      -> DrawNormalized();
    Response_UP              -> DrawNormalized("same");
    Response_DOWN        -> DrawNormalized("same");
    Response_STRANGE   -> DrawNormalized("same");
    Response_CHARM      -> DrawNormalized("same");
    Response_BOTTOM    -> DrawNormalized("same");
    leg2 -> Draw();
    canvas ->SaveAs("Histogram_SmearingFunctions/BinPt200GeV/SixCategory/Compare_Response.png");
     
  }//end if bool

  /////////////////////////////////////////////////////////////////////////////////

  if( Draw_SmearingFunctions_Bin200_FourCategory == true){
     
    TFile file1("SmearingFunctions_RecoGen_GravitonSamples/BinPt200GeV/FourCategory/SmearingFunctions_ResonanceSamplesToQuarkQuark_Bin200.root ");
    TH1D *H_step_pt = (TH1D*)file1.Get("H_step_pt");
    int step_pt = H_step_pt->GetMean();
    TFile file2("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/FourCategory/SmearingFunctions_QCDSamples_Bin200.root ");    
    
    int n_bin = 4500 / step_pt ;   
    char HistoName[200];	 
    TH1D *h[100][100];
    TH1D *j[100][100];
    
    TLegend *leg = new TLegend( 0.9, 0.9, 0.6, 0.6);
    leg -> SetHeader("Smearing Functions: ");
    leg -> AddEntry(h[0][0], "Resonance Samples", "L");
    leg -> AddEntry(j[0][0], "QCD Samples", "L");
    
    for( int kk=0; kk<4; kk++){//tipi di partone   
      for(int ii=0; ii<5; ii++){//bin in eta
	for(int jj=0; jj<n_bin ; jj++){//bin in pT	 	 
	  int pt_bin_max = (jj*step_pt)+step_pt;
	  // da cambiare e metterci i nomi dei quark
	  if( kk == 0 )	sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 1 )	sprintf(HistoName,"Histo_RQuarkLIGHT_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 2 )	sprintf(HistoName,"Histo_RQuarkCHARM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  if( kk == 3 )	sprintf(HistoName,"Histo_RQuarkBOTTOM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  j[ii][jj] = (TH1D*)file2.Get(HistoName);   
	  cout<<"      Graviton->Integral(): "<<       h[ii][jj]->Integral() <<endl;
	  cout<<"      QCD->Integral(): "<<       j[ii][jj]->Integral() <<endl;
	  Normalizer(h[ii][jj]);
	  Normalizer( j[ii][jj]);
	  char output_Dir[200] = "Histogram_SmearingFunctions/BinPt200GeV/FourCategory/";
	  char extension[10] = ".png";
	  strcat( HistoName, extension);
	  cout<<HistoName<<endl;
	  if( j[ii][jj]->Integral() >= h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, j[ii][jj], h[ii][jj], 0., 2.,"R","Normalized", kBlue, kRed);
	  if( j[ii][jj]->Integral() < h[ii][jj]->Integral() )	DrawAndSave(output_Dir, HistoName, h[ii][jj], j[ii][jj], 0., 2.,"R","Normalized", kRed, kBlue );

	}
      }
    }

    TH1D *Response_LIGHT       = (TH1D*)file2.Get("Histo_RQuarkLIGHT_WideJet_GenWideJet_0_1600");   
    TH1D *Response_CHARM    = (TH1D*)file2.Get("Histo_RQuarkCHARM_WideJet_GenWideJet_0_1600");   
    TH1D *Response_BOTTOM  = (TH1D*)file2.Get("Histo_RQuarkBOTTOM_WideJet_GenWideJet_0_1600");   
    TH1D *Response_GLUON    = (TH1D*)file2.Get("Histo_RGluons_WideJet_GenWideJet_0_1600");   
    
    TCanvas *canvas = new TCanvas("canvas","",800,800);
    Response_GLUON    -> SetTitle(" ");
    Response_GLUON    -> SetStats(0);
    Response_GLUON     ->GetXaxis()->SetRangeUser(0.4, 1.3);
    Response_GLUON    -> SetXTitle("Response");
    Response_GLUON    -> SetYTitle("Normalized to unit");
    Response_GLUON    -> GetYaxis()->SetTitleOffset(1.55);
    Response_LIGHT         -> SetLineColor(kRed);
    Response_CHARM       -> SetLineColor(kViolet);
    Response_BOTTOM     -> SetLineColor(kYellow+1);
    Response_GLUON       -> SetLineColor(kBlack);    
    Response_LIGHT         -> SetMarkerColor(kRed);
    Response_CHARM       -> SetMarkerColor(kViolet);
    Response_BOTTOM     -> SetMarkerColor(kYellow+1);
    Response_GLUON       -> SetMarkerColor(kBlack);    
    Response_LIGHT         -> SetMarkerStyle(20);
    Response_CHARM      -> SetMarkerStyle(23);
    Response_BOTTOM    -> SetMarkerStyle(31);
    Response_GLUON       -> SetMarkerStyle(33);   
    TLegend *leg2 = new TLegend( 0.1, 0.9, 0.4, 0.6);
    leg2 -> SetHeader("Parton Flavour:");
    leg2 -> AddEntry(Response_LIGHT, "quark Light", "L");
    leg2 -> AddEntry(Response_CHARM, "quark Charm", "L");
    leg2 -> AddEntry(Response_BOTTOM, "quark Bottom", "L");
    leg2 -> AddEntry(Response_GLUON, "Gluons", "L");
    Response_GLUON      -> DrawNormalized();
    Response_LIGHT        -> DrawNormalized("same");
    Response_CHARM      -> DrawNormalized("same");
    Response_BOTTOM    -> DrawNormalized("same");
    leg2 -> Draw();
    canvas ->SaveAs("Histogram_SmearingFunctions/BinPt200GeV/FourCategory/Compare_Response.png");

  }//end if bool

  //////////////////////////////////////////////////////////////////////////////////////////////

  bool Modifica = false;

  if(Modifica == true){
    TFile file5("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/SixCategory/SmearingFunctions_QCDSamples_Bin200.root ");
    TFile file6("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/FourCategory/SmearingFunctions_QCDSamples_Bin200.root ");    

    TH1D *H_step_pt = (TH1D*)file5.Get("H_step_pt");
    int step_pt = H_step_pt->GetMean();

    
    int n_bin = 4500 / step_pt ;   
    char HistoName[200];	 
    TH1D *h[100][100];
    TH1D *j[100][100];
    
    TLegend *leg = new TLegend( 0.9, 0.9, 0.6, 0.6);
    leg -> SetHeader("Smearing Functions: ");
    leg -> AddEntry(h[0][0], "Six Categories", "L");
    leg -> AddEntry(j[0][0], "Four Categories", "L");

      for(int ii=0; ii<5; ii++){//bin in eta
	for(int jj=0; jj<n_bin ; jj++){//bin in pT	 	 
	  int pt_bin_max = (jj*step_pt)+step_pt;
	  sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  h[ii][jj]= (TH1D*)file5.Get(HistoName);   
	  j[ii][jj] = (TH1D*)file6.Get(HistoName);   
	  cout<<"      Six categories->Integral(): "<<       h[ii][jj]->Integral() <<endl;
	  cout<<"      Four categories->Integral(): "<<       j[ii][jj]->Integral() <<endl;

	}
      }
  }//end if bool

  ////////////////////////////////////////////////////////////////////////////////

  if( Draw_SmearingFunctions_Bin200_FourCategory_Alone == true){
     
  TFile file1("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/FourCategory/SmearingFunctions_QCDSamples_Bin200.root ");    
  TH1D *H_step_pt = (TH1D*)file1.Get("H_step_pt");
  int step_pt = H_step_pt->GetMean();
     
  int n_bin = 4500 / step_pt ;   
  char HistoName[200];	 
  TH1D *h[100][100];
  TH1D *j[100][100];
      
  //    for( int kk=0; kk<4; kk++){//tipi di partone   
      for(int ii=0; ii<5; ii++){//bin in eta
	for(int jj=0; jj<n_bin ; jj++){//bin in pT	 	 
	  int pt_bin_max = (jj*step_pt)+step_pt;




	  //	  cout<<" "<< endl;
	  sprintf(HistoName,"Histo_RGluons_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  //	  cout<<HistoName<<"        ->Mean(): "<<       h[ii][jj]->GetMean() <<" +/- "<<h[ii][jj]->GetRMS() <<endl;




	  sprintf(HistoName,"Histo_RQuarkLIGHT_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  //	  cout<<HistoName<<"->Mean(): "<<       h[ii][jj]->GetMean() <<" +/- "<<h[ii][jj]->GetRMS() <<endl;





	  // sprintf(HistoName,"Histo_RQuarkCHARM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  // h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  // cout<<HistoName<<"->Mean(): "<<       h[ii][jj]->GetMean() <<" +/- "<<h[ii][jj]->GetRMS() <<endl;
	  // sprintf(HistoName,"Histo_RQuarkBOTTOM_WideJet_GenWideJet_%d_%d",ii,pt_bin_max);
	  // h[ii][jj]= (TH1D*)file1.Get(HistoName);   
	  // cout<<HistoName<<"->Mean(): "<<       h[ii][jj]->GetMean() <<" +/- "<<h[ii][jj]->GetRMS() <<endl;
	  

	}
      }

  }//end if bool
  ////////////////////////////////////////////////////////////////////////////



  if( Draw_SmearingFunctions_Bin200_TwoCategory_Example == true){
     
  TFile file1("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/TwoCategory/SmearingFunctions_QCDSamples_Bin200.root ");    
  TH1D *Histo_RQuarks_WideJet_GenWideJet_0_1600 = (TH1D*)file1.Get("Histo_RQuarks_WideJet_GenWideJet_0_1600");
  TH1D *Histo_RQuarks_WideJet_GenWideJet_2_1200 = (TH1D*)file1.Get("Histo_RQuarks_WideJet_GenWideJet_2_1200");

  char output_Dir[200] = "Histogram_SmearingFunctions/BinPt200GeV/TwoCategory/";
  char extension[10] = ".png";
  char HistoName[200];	 

  //voglio fare fit gaussiani a questi due

  TF1 *fgaus = new TF1("fgaus","gaus", 0.9, 1.1);
  //  fgaus->SetParLimits(0,0,300000);
  //fgaus->SetParLimits(1,-0.2,0.2);   
  //fgaus->SetParLimits(2,-0.1,0.1);

  double N, Mean, Sigma;
  double min_fit, max_fit;

    N        = Histo_RQuarks_WideJet_GenWideJet_0_1600->GetEntries();        
    Mean  = Histo_RQuarks_WideJet_GenWideJet_0_1600->GetMean();           
    Sigma = Histo_RQuarks_WideJet_GenWideJet_0_1600->GetRMS();            

    // min_fit = Mean - 2*Sigma;
    // max_fit = Mean + 2*Sigma;

     min_fit = 0.93;
     max_fit = 1.09;

    fgaus ->SetParameters(N, Mean, Sigma);   
    Histo_RQuarks_WideJet_GenWideJet_0_1600->Fit("fgaus","lR", " ", min_fit, max_fit);
    
    sprintf(HistoName,"Histo_RQuarks_WideJet_GenWideJet_0_1600");
    strcat( HistoName, extension);
    DrawAndSave(output_Dir, HistoName, Histo_RQuarks_WideJet_GenWideJet_0_1600, 0.7 , 1.3,"R","Events");


    N        = Histo_RQuarks_WideJet_GenWideJet_2_1200->GetEntries();        
    Mean  = Histo_RQuarks_WideJet_GenWideJet_2_1200->GetMean();           
    Sigma = Histo_RQuarks_WideJet_GenWideJet_2_1200->GetRMS();            

     min_fit = 0.91;
     max_fit = 1.13;

    fgaus ->SetParameters(N, Mean, Sigma);   
    Histo_RQuarks_WideJet_GenWideJet_2_1200->Fit("fgaus","lR", " ", min_fit, max_fit);

    sprintf(HistoName,"Histo_RQuarks_WideJet_GenWideJet_2_1200");
    strcat( HistoName, extension);
    DrawAndSave(output_Dir, HistoName, Histo_RQuarks_WideJet_GenWideJet_2_1200, 0.7 , 1.3,"R","Events");

  //  DrawAndSave(output_Dir, "Histo_RQuarks_WideJet_GenWideJet_2_1200", Histo_RQuarks_WideJet_GenWideJet_2_1200, 0., 2.,"R","Events");


}//end if bool


  TFile file1("SmearingFunctions_RecoGen_QCDSamples/BinPt200GeV/FourCategory/SmearingFunctions_QCDSamples_Bin200.root ");    
  TH1D *HistoQCD_GenWideJet1_Pt = (TH1D*)file1.Get("H_GenWideJet1_pT");
  TH1D *HistoQCD_GenWideJet1_Eta = (TH1D*)file1.Get("H_GenWideJet1_Eta");

  TFile file2("SmearingFunctions_RecoGen_GravitonSamples/BinPt200GeV/FourCategory/SmearingFunctions_ResonanceSamplesToQuarkQuark_Bin200.root");
  TH1D *HistoRSG_GenWideJet1_Pt = (TH1D*)file2.Get("H_GenWideJet1_pT");

  Normalizer(HistoRSG_GenWideJet1_Pt);
  Normalizer(HistoQCD_GenWideJet1_Pt);

  TLegend *leg = new TLegend( 0.9, 0.9, 0.6, 0.6);
  leg -> SetHeader(" ");
  leg -> AddEntry(HistoQCD_GenWideJet1_Pt, "QCD Samples", "L");
  leg -> AddEntry(HistoRSG_GenWideJet1_Pt, "Signal Samples", "L");


  TCanvas *canvas1 = new TCanvas("canvas1","",800,800);
  HistoQCD_GenWideJet1_Pt    -> SetTitle(" ");
  HistoQCD_GenWideJet1_Pt    -> SetLineColor(kBlue);
  HistoQCD_GenWideJet1_Pt    -> SetStats(0);
  HistoQCD_GenWideJet1_Pt    -> SetXTitle("Pt [GeV]");
  HistoQCD_GenWideJet1_Pt    -> SetYTitle("Events");
  HistoQCD_GenWideJet1_Pt    -> GetYaxis()->SetTitleOffset(1.45);
  canvas1->SetLogy();
  HistoQCD_GenWideJet1_Pt    -> GetXaxis()->SetRangeUser(0, 5000);
  HistoQCD_GenWideJet1_Pt    -> Draw();
  canvas1 ->SaveAs("QCDGenWide_Pt.png");



  TCanvas *canvas2 = new TCanvas("canvas2","",800,800);
  HistoQCD_GenWideJet1_Eta    -> SetTitle(" ");
  HistoQCD_GenWideJet1_Eta    -> SetLineColor(kBlue);
  HistoQCD_GenWideJet1_Eta    -> SetStats(0);
  HistoQCD_GenWideJet1_Eta    -> SetXTitle("#eta");
  HistoQCD_GenWideJet1_Eta    -> SetYTitle("Events");
  HistoQCD_GenWideJet1_Eta    -> GetYaxis()->SetTitleOffset(1.45);
  HistoQCD_GenWideJet1_Eta    -> GetXaxis()->SetRangeUser(-2.6, 2.6);
  HistoQCD_GenWideJet1_Eta    -> Draw();
  canvas2 ->SaveAs("QCDGenWide_Eta.png");


  ///////////////////////////////////////////////////////////////////////////

}//end main
