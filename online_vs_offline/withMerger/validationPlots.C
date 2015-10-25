#include "TH1D.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "TColor.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAttPad.h"
#include <iostream>

void validationPlots()
{
  bool isMB = 1;  

  int doCentCut = 0;
  int doEtaCut = 0; 
  int doPtCut = 0;
  const char * ptCut = "trkPt>30 && trkPt<300";
  int doFake = 0;
  int doVtx = 1;
  int nEvt = 2;

  bool doFineMVABins = 0;

  TH1::SetDefaultSumw2();
  TLatex * lat = new TLatex(0.5,0.5,"test");

  TFile * fMC;
  TFile * fDa;
  //if(isMB) fMC = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/hiForest_HydjetMB_2076GeV_FOREST_753p1_merged/HydjetMB_2076GeV_FOREST_753p1_v0_merged.root","read");
  if(isMB) fDa = TFile::Open("HiForest.root","read");
  if(!isMB) fMC = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/mergedForest/HiForest_HydjetMB_Pyquen_DiJet_pt80to9999_2670GeV_cfi_753_patch1/HiForest_HydjetMB_Pyquen_DiJet_pt80to9999_2670GeV_cfi_753_patch1_run1_mc_v2.root","read");
  if(!isMB) fDa = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/HiForest_HIHighPt_HIRun2011-v1_RECO_753_patch1_Jet80_merged/HiForest_HIHighPt_jet80_HIRun2011-v1_merged.root","read");

  TTree * tree[2];
  tree[0] = (TTree*)fDa->Get("hltAnaTrack/trackTree");
  tree[1] = (TTree*)fDa->Get("anaTrack/trackTree");
  TTree * hiBinTree[2];
  hiBinTree[0] = (TTree*)fDa->Get("hiEvtAnalyzer/HiTree");
  hiBinTree[1] = (TTree*)fDa->Get("hiEvtAnalyzer/HiTree");

  TH1D * daVtx = new TH1D("daVtx",";z",30,-15,15);
  TH1D * mcVtx = new TH1D("mcVtx",";z",30,-15,15);
  tree[1]->Draw("zVtx[0]>>daVtx","","",100000);
  tree[0]->Draw("zVtx[0]>>mcVtx","","",100000);
  daVtx->Scale(1.0/daVtx->Integral(1,30));
  mcVtx->Scale(1.0/mcVtx->Integral(1,30));
 // mcVtx->Divide(daVtx);

  TH1D * daHiBin = new TH1D("daHiBin",";hiBin",100,0,200);
  TH1D * mcHiBin = new TH1D("mcHiBin",";hiBin",100,0,200);
  hiBinTree[1]->Draw("hiBin>>daHiBin","");
  hiBinTree[0]->Draw("hiBin>>mcHiBin","");
  daHiBin->Scale(1.0/daHiBin->Integral(1,100));
  mcHiBin->Scale(1.0/mcHiBin->Integral(1,100));
//  mcHiBin->Divide(daHiBin);
  std::cout << "here" << std::endl;


  int hiBin = 0;
  float hiBinw = 0;
  float vz[200] = {0};
  float weight = 0;
  tree[1]->SetBranchAddress("zVtx",&vz);
  hiBinTree[1]->SetBranchAddress("hiBin",&hiBin);
  TFile * f = TFile::Open("vertexWeight.root","recreate");
  TTree * w = new TTree("w","w");
  w->Branch("vtxw",&weight);
  w->Branch("hiBinw",&hiBinw);
  for(int i=0; i<nEvt; i++)
  {
    if(i%1000==0) std::cout << i << std::endl;
    tree[1]->GetEntry(i);
    hiBinTree[1]->GetEntry(i);
    weight=mcVtx->GetBinContent(mcVtx->FindBin(vz[0]));
    hiBinw=mcHiBin->GetBinContent(mcHiBin->FindBin(hiBin));
    weight=1;
    hiBinw=1;
    w->Fill();
  }
  w->Write();
  tree[1]->AddFriend(w);
  tree[1]->AddFriend(hiBinTree[1]);
  tree[0]->AddFriend(w);
  tree[0]->AddFriend(hiBinTree[0]);
  std::cout << "done writing vtx tree" << std::endl;

  TH1D *chi2[6][2][2], *dxy[6][2][2], *dz[6][2][2], *nhit[6][2][2], *nlayer[6][2][2], *eta[6][2][2], *pterr[6][2][2], *mva[6][2][2], *mvaRat[6][2][1];

  int algos[5] = {4,5,6,7,11};
  const char *algName[5] = {" (Initial)"," (low pt triplet)"," (pixel pair)"," (detached)"," (jet-core)"};
  int nTracks[6][2][2];

  for(int algo = 0; algo<6; algo++)
  {
    for(int purity = 0; purity <2; purity++)
    {
      for(int sample = 0; sample<2; sample++)
      {
        chi2[algo][purity][sample] = new TH1D(Form("chi2%d%d%d",algo,purity,sample),";chi2/ndof/nlayer;dN/dchi2",100,0,1);
        chi2[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkChi2/(1.0*trkNlayer*trkNdof)>>chi2%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d)) && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d)) && ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        chi2[algo][purity][sample]->Scale(1.0/chi2[algo][purity][sample]->Integral(1,100));
        chi2[algo][purity][sample]->GetYaxis()->SetRangeUser(0,0.2);
        std::cout << algo << std::endl;  
    
        dxy[algo][purity][sample] = new TH1D(Form("dxy%d%d%d",algo,purity,sample),";dxy/dxyerr;dN/(dxy/dxyerr)",100,-100,100);
        dxy[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkDxy1/(1.0*trkDxyError1)>>dxy%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        dxy[algo][purity][sample]->Scale(1.0/dxy[algo][purity][sample]->Integral(1,100));
        dxy[algo][purity][sample]->GetYaxis()->SetRangeUser(10e-7,1);
        
        dz[algo][purity][sample] = new TH1D(Form("dz%d%d%d",algo,purity,sample),";dz/dzerr;dN/(dz/dzerr)",100,-100,100);
        dz[algo][purity][sample]->SetMarkerSize(0.8);
        tree[sample]->Draw(Form("trkDz1/(1.0*trkDzError1)>>dz%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        dz[algo][purity][sample]->Scale(1.0/dz[algo][purity][sample]->Integral(1,100));
        dz[algo][purity][sample]->GetYaxis()->SetRangeUser(10e-7,1);

        nhit[algo][purity][sample] = new TH1D(Form("nhit%d%d%d",algo,purity,sample),";nhit;dN/nhit",28,3,30);
        nhit[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkNHit>>nhit%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        nhit[algo][purity][sample]->Scale(1.0/nhit[algo][purity][sample]->Integral(1,28));
        nhit[algo][purity][sample]->GetYaxis()->SetRangeUser(0,0.25);
        
        nlayer[algo][purity][sample] = new TH1D(Form("nlayer%d%d%d",algo,purity,sample),";nlayer;dN/nlayer",22,4,25);
        nlayer[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkNlayer>>nlayer%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        nlayer[algo][purity][sample]->Scale(1.0/nlayer[algo][purity][sample]->Integral(1,22));
        nlayer[algo][purity][sample]->GetYaxis()->SetRangeUser(0.0,0.3);
        
        eta[algo][purity][sample] = new TH1D(Form("eta%d%d%d",algo,purity,sample),";eta;dN/eta",50,-2.4,2.4);
        eta[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkEta>>eta%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        eta[algo][purity][sample]->Scale(1.0/eta[algo][purity][sample]->Integral(1,50));
        eta[algo][purity][sample]->GetYaxis()->SetRangeUser(0.0,0.03);
        
        pterr[algo][purity][sample] = new TH1D(Form("pterr%d%d%d",algo,purity,sample),";pterr;dN/(pt/pterr)",50,0,0.05);
        pterr[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkPtError/trkPt>>pterr%d%d%d",algo,purity,sample),Form("(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        pterr[algo][purity][sample]->Scale(1.0/pterr[algo][purity][sample]->Integral(1,50));
        pterr[algo][purity][sample]->GetYaxis()->SetRangeUser(0.0,0.15);
        std::cout << algo << std::endl;  
        
        mva[algo][purity][sample] = new TH1D(Form("mva%d%d%d",algo,purity,sample),";mva;dN/d(mva)",doFineMVABins?100:30,-1,1);
        mva[algo][purity][sample]->SetMarkerSize(0.8); 
        tree[sample]->Draw(Form("trkMVA>>mva%d%d%d",algo,purity,sample),Form("(1+0.2*(!(%d) && %d && trkFake))*(1+(vtxw-1)*(%d&&%d))*(1+(hiBinw-1)*(%d))*(((highPurity== 1)||(highPurity==%d))  && (trkEta>-2.4 && trkEta<2.4) && (%d==5 || %d==trkAlgo) && ((%s) || (!%d))&& ((trkEta>0.8 || trkEta<-0.8) || (!%d))&& ((hiBin>=100 && hiBin<200) || (!%d)))",sample,doFake,sample,doVtx,sample,purity,algo,algos[algo],ptCut,doPtCut,doEtaCut,doCentCut),"",nEvt);
        nTracks[algo][purity][sample]=mva[algo][purity][sample]->Integral(1,doFineMVABins?100:30);
        mva[algo][purity][sample]->Scale(1.0/mva[algo][purity][sample]->Integral(1,doFineMVABins?100:30));
        mva[algo][purity][sample]->GetYaxis()->SetRangeUser(0.0,0.4);
        std::cout << algo << std::endl;  

        if(sample==1)
        {
          mvaRat[algo][purity][0] = (TH1D*)mva[algo][purity][1]->Clone(Form("mvaRat%d%d0",algo,purity));
          mvaRat[algo][purity][0]->GetYaxis()->SetTitle("offline/online");
          mvaRat[algo][purity][0]->Divide(mva[algo][purity][0]);
          mvaRat[algo][purity][0]->GetYaxis()->SetRangeUser(0,2);
          std::cout << algo << std::endl;  
        }
      }
    }
  }

  TCanvas * c1[2];
  TLegend * l1[2];
  for(int j = 0; j<2; j++)
  {
    c1[j] = new TCanvas(Form("c1_%d",j),Form("c1_%d",j),1800,1000);
    c1[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue;
      c1[j]->cd(i+1);
      chi2[i][j][0]->Draw("p");

      chi2[i][j][1]->SetMarkerColor(kRed+1);
      chi2[i][j][1]->SetLineColor(kRed+1);
      chi2[i][j][1]->Draw("p same");

      if(i==0)
      {
        l1[j] = new TLegend(0.6,0.5,0.9,0.9);
        l1[j]->AddEntry(chi2[0][j][0],"Online");
        l1[j]->AddEntry(chi2[0][j][1],"Offline");
        if(j==1) l1[j]->AddEntry((TObject*)0,"HighPurity","");
        l1[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(0.5,0.04,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(0.5,0.04,"All Algos");
    } 
    c1[j]->SaveAs(Form("validationPlots/chi2_%d.png",j));
  }
  
  TCanvas * c2[2];
  TLegend * l2[2];
  for(int j = 0; j<2; j++)
  {
    c2[j] = new TCanvas(Form("c2_%d",j),Form("c2_%d",j),1800,1000);
    c2[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue;
      c2[j]->cd(i+1);
      c2[j]->cd(i+1)->SetLogy();
      dxy[i][j][0]->Draw("p");

      dxy[i][j][1]->SetMarkerColor(kRed+1);
      dxy[i][j][1]->SetLineColor(kRed+1);
      dxy[i][j][1]->Draw("p same");

      if(i==0)
      {
        l2[j] = new TLegend(0.6,0.5,0.9,0.9);
        l2[j]->AddEntry(dxy[0][j][0],"Online");
        l2[j]->AddEntry(dxy[0][j][1],"Offline");
        if(j==1) l2[j]->AddEntry((TObject*)0,"HighPurity","");
        l2[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(-80,0.1,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(-80,0.1,"All Algos");
    }
    c2[j]->SaveAs(Form("validationPlots/dxy_%d.png",j));
  }
  
  TCanvas * c3[2];
  TLegend * l3[2];
  for(int j = 0; j<2; j++)
  {
    c3[j] = new TCanvas(Form("c3_%d",j),Form("c3_%d",j),1800,1000);
    c3[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue;
      c3[j]->cd(i+1);
      c3[j]->cd(i+1)->SetLogy();
      dz[i][j][0]->Draw("p");

      dz[i][j][1]->SetMarkerColor(kRed+1);
      dz[i][j][1]->SetLineColor(kRed+1);
      dz[i][j][1]->Draw("p same");

      if(i==0)
      {
        l3[j] = new TLegend(0.6,0.5,0.9,0.9);
        l3[j]->AddEntry(dz[0][j][0],"Online");
        l3[j]->AddEntry(dz[0][j][1],"Offline");
        if(j==1) l3[j]->AddEntry((TObject*)0,"HighPurity","");
        l3[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(-80,0.1,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(-80,0.1,"All Algos");
    }
    c3[j]->SaveAs(Form("validationPlots/dz_%d.png",j));
  }
  
  TCanvas * c4[2];
  TLegend * l4[2];
  for(int j = 0; j<2; j++)
  {
    c4[j] = new TCanvas(Form("c4_%d",j),Form("c4_%d",j),1800,1000);
    c4[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue;
      c4[j]->cd(i+1);
      nhit[i][j][0]->Draw("p");

      nhit[i][j][1]->SetMarkerColor(kRed+1);
      nhit[i][j][1]->SetLineColor(kRed+1);
      nhit[i][j][1]->Draw("p same");

      if(i==0)
      {
        l4[j] = new TLegend(0.6,0.5,0.9,0.9);
        l4[j]->AddEntry(nhit[0][j][0],"Online");
        l4[j]->AddEntry(nhit[0][j][1],"Offline");
        if(j==1) l4[j]->AddEntry((TObject*)0,"HighPurity","");
        l4[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(5,0.15,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(5,0.15,"All Algos");
    }
    c4[j]->SaveAs(Form("validationPlots/nhit_%d.png",j));
  }
  
  TCanvas * c5[2];
  TLegend * l5[2];
  for(int j = 0; j<2; j++)
  {
    c5[j] = new TCanvas(Form("c5_%d",j),Form("c5_%d",j),1800,1000);
    c5[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue; 
      c5[j]->cd(i+1);
      nlayer[i][j][0]->Draw("p");

      nlayer[i][j][1]->SetMarkerColor(kRed+1);
      nlayer[i][j][1]->SetLineColor(kRed+1);
      nlayer[i][j][1]->Draw("p same");

      if(i==0)
      {
        l5[j] = new TLegend(0.6,0.5,0.9,0.9);
        l5[j]->AddEntry(nlayer[0][j][0],"Online");
        l5[j]->AddEntry(nlayer[0][j][1],"Offline");
        if(j==1) l5[j]->AddEntry((TObject*)0,"HighPurity","");
        l5[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(5,0.2,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(5,0.2,"All Algos");
    }
    c5[j]->SaveAs(Form("validationPlots/nlayer_%d.png",j));
  }
  
  TCanvas * c6[2];
  TLegend * l6[2];
  for(int j = 0; j<2; j++)
  {
    c6[j] = new TCanvas(Form("c6_%d",j),Form("c6_%d",j),1800,1000);
    c6[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue; 
      c6[j]->cd(i+1);
      eta[i][j][0]->Draw("p");

      eta[i][j][1]->SetMarkerColor(kRed+1);
      eta[i][j][1]->SetLineColor(kRed+1);
      eta[i][j][1]->Draw("p same");

      if(i==0)
      {
        l6[j] = new TLegend(0.5,0.2,0.9,0.5);
        l6[j]->AddEntry(eta[0][j][0],"Online");
        l6[j]->AddEntry(eta[0][j][1],"Offline");
        if(j==1) l6[j]->AddEntry((TObject*)0,"HighPurity","");
        l6[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(-1.5,0.005,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(-1.5,0.005,"All Algos");
    }
    c6[j]->SaveAs(Form("validationPlots/eta_%d.png",j));
  }
  
  TCanvas * c7[2];
  TLegend * l7[2];
  for(int j = 0; j<2; j++)
  {
    c7[j] = new TCanvas(Form("c7_%d",j),Form("c7_%d",j),1800,1000);
    c7[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue;
      c7[j]->cd(i+1);
      pterr[i][j][0]->Draw("p");

      pterr[i][j][1]->SetMarkerColor(kRed+1);
      pterr[i][j][1]->SetLineColor(kRed+1);
      pterr[i][j][1]->Draw("p same");

      if(i==0)
      {
        l7[j] = new TLegend(0.6,0.5,0.9,0.9);
        l7[j]->AddEntry(pterr[0][j][0],"Online");
        l7[j]->AddEntry(pterr[0][j][1],"Offline");
        if(j==1) l7[j]->AddEntry((TObject*)0,"HighPurity","");
        l7[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(0.03,0.04,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(0.03,0.04,"All Algos");
    }
    c7[j]->SaveAs(Form("validationPlots/pterr_%d.png",j));
  }
 
  TCanvas * c8[2];
  TLegend * l8[2];
  for(int j = 0; j<2; j++)
  {
    c8[j] = new TCanvas(Form("c8_%d",j),Form("c8_%d",j),1800,1000);
    c8[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4 && isMB) continue;
      c8[j]->cd(i+1);
      mva[i][j][0]->Draw("p");

      mva[i][j][1]->SetMarkerColor(kRed+1);
      mva[i][j][1]->SetLineColor(kRed+1);
      mva[i][j][1]->Draw("p same");

      if(i==0)
      {
        l8[j] = new TLegend(0.6,0.5,0.9,0.9);
        l8[j]->AddEntry(mva[0][j][0],"Online");
        l8[j]->AddEntry(mva[0][j][1],"Offline");
        if(j==1) l8[j]->AddEntry((TObject*)0,"HighPurity","");
        l8[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(0.03,0.04,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(0.03,0.04,"All Algos");
    }
    c8[j]->SaveAs(Form("validationPlots/mva_%d.png",j));
  }
  
  TCanvas * c9[2];
  TLegend * l9[2];
  for(int j = 0; j<2; j++)
  {
    c9[j] = new TCanvas(Form("c9_%d",j),Form("c9_%d",j),1800,1000);
    c9[j]->Divide(3,2);
    for(int i = 0; i<6; i++)
    {
      if(i==4) continue;// && isMB) continue;
      c9[j]->cd(i+1);
      mvaRat[i][j][0]->Draw("p");

    //  mvaRat[i][j][1]->SetMarkerColor(kRed+1);
    //  mvaRat[i][j][1]->SetLineColor(kRed+1);
    //  mvaRat[i][j][1]->Draw("p same");

      if(i==0)
      {
        l9[j] = new TLegend(0.6,0.5,0.9,0.9);
      //  l9[j]->AddEntry(mvaRat[0][j][0],"Online");
      //  l9[j]->AddEntry(mvaRat[0][j][1],"Offline");
        if(j==1) l9[j]->AddEntry((TObject*)0,"HighPurity","");
        l9[j]->Draw("same");
      }
      if(i<5) lat->DrawLatex(0.03,0.4,Form("Algo %d%s",algos[i],algName[i]));
      if(i==5) lat->DrawLatex(0.03,0.4,"All Algos");
    }
    c9[j]->SaveAs(Form("validationPlots/mvaRat_%d.png",j));
  }

  //change in signal by reweighting MVAs
  for(int i = 0; i<6; i++)
  {
    if(i==4) continue;
    float cut[6] = {-0.77,0.35,0.77,-0.09,0,0};
    float dSig = mva[i][0][1]->Integral(mva[i][0][1]->FindBin(cut[i]),mva[i][0][1]->FindBin(1)); 
    mva[i][0][1]->Divide(mvaRat[i][0][0]);
    float dSig2 = (float)mva[i][0][1]->Integral(mva[i][0][1]->FindBin(cut[i]),mva[i][0][1]->FindBin(1))/mva[i][0][1]->Integral(mva[i][0][1]->FindBin(-1),mva[i][0][1]->FindBin(1)); 
    std::cout << "\n\nAlgo " << algos[i] << "\nSignal Purity before MVA weight: " << dSig << "  Signal Purity after MVA weight: " << dSig2 << "  Difference (as %): " << dSig-dSig2 << " (" << (dSig-dSig2)*100/dSig<<")" << std::endl;
    std::cout << "Fraction of tracks w/o highPurity (data): " << (float)nTracks[i][0][1]/nTracks[5][0][1] << std::endl;
    std::cout << "Fraction of tracks w/o highPurity (MC): " << (float)nTracks[i][0][0]/nTracks[5][0][0] << std::endl;
    std::cout << "\nFraction of tracks w/  highPurity  (data): " << (float)nTracks[i][1][1]/nTracks[5][1][1] << std::endl;
    std::cout << "Fraction of tracks w/  highPurity (MC): " << (float)nTracks[i][1][0]/nTracks[5][1][0] << std::endl;
  }
}
