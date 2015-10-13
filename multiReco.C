#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TAttLine.h"
#include <iostream>

void multiReco()
{
  TH1::SetDefaultSumw2();
  TH1D * num50 = new TH1D("num50",";pPt;multiReco",30,0.5,200);
  TH1D * den50 = new TH1D("den50",";pPt;multiReco",30,0.5,200);
  TFile * f = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/hiForest_HydjetNcoll_Pyquen_DiJet_pt220to9999_5020GeV_cfi_750_run2_mc_HIon_50percent_merged/hiForest_HydjetNcoll_Pyquen_DiJet_pt220to9999_5020GeV_cfi_750_run2_mc_HIon_50percent_multrec.root","Read"); 
  TTree * m50 = (TTree*)f->Get("anaTrack/trackTree");

  int nPart; 
  int pNRec[100000];
  float pt[100000];

  m50->SetBranchAddress("nParticle",&nPart);
  m50->SetBranchAddress("pPt", &pt);
  m50->SetBranchAddress("pNRec", &pNRec);

  for(int i = 0; i<1000; i++)//i<m50->GetEntries();i++)
  {
    if(i%1000 == 0) std::cout << i << "/" << m50->GetEntries() << std::endl;
    m50->GetEntry(i);
    for(int j = 0; j<nPart; j++)
    {
      den50->Fill(pt[j]);
      if(pNRec[j]>1) num50->Fill(pt[j],pNRec[j]-1); 
    }
  }
  den50->Print("All"); num50->Print("All");
  num50->Divide(den50);
  num50->Print("All");
  f->Close();

 //*************************************************************************************
  TH1D * num75 = new TH1D("num75",";pPt;multiReco",30,0.5,200);
  TH1D * den75 = new TH1D("den75",";pPt;multiReco",30,0.5,200);
  f = TFile::Open("/mnt/hadoop/cms/store/user/dgulhan/hiForest_HydjetNcoll_Pyquen_DiJet_pt220to9999_5020GeV_cfi_750_run2_mc_HIon_BS_calomatching/hiForest_HydjetNcoll_Pyquen_DiJet_pt220to9999_5020GeV_cfi_750_run2_mc_HIon_BS_calomatching_merged_forest_3_0.root","Read"); 
  m50 = (TTree*)f->Get("anaTrack/trackTree");

  m50->SetBranchAddress("nParticle",&nPart);
  m50->SetBranchAddress("pPt", &pt);
  m50->SetBranchAddress("pNRec", &pNRec);

  for(int i = 0;i<1000;i++)// i<m50->GetEntries();i++)
  {
    if(i%1000 == 0) std::cout << i << "/" << m50->GetEntries() << std::endl;
    m50->GetEntry(i);
    for(int j = 0; j<nPart; j++)
    {
      den75->Fill(pt[j]);
      if(pNRec[j]>1) num75->Fill(pt[j],pNRec[j]-1); 
    }
  }
  den75->Print("All"); num75->Print("All");
  num75->Divide(den75);
  num75->Print("All");

  num50->Draw("p");
  num75->SetLineColor(kRed);
  num75->Draw("h same");
}
