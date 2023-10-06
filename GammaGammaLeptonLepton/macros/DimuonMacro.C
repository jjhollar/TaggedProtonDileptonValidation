#define DimuonMacro_cxx
#include "DimuonMacro.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void DimuonMacro::Loop(Int_t multi, Int_t mc, Int_t sb, Int_t yr, Int_t nearfar)
{
//   In a ROOT session, you can do:
//      root> .L DimuonMacro.C
//      root> DimuonMacro t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   TH1F *hacop = new TH1F("hacop","hacop",1000,0,0.01);
   TH2F *hcorr45 = new TH2F("hcorr45","hcorr45",500,0,0.25,500,0,0.25);
   TH2F *hcorr56 = new TH2F("hcorr56","hcorr56",500,0,0.25,500,0,0.25);

   TH2F *hcorrmult45 = new TH2F("hcorrmult45","hcorrmult45",500,0,0.25,500,0,0.25);
   TH2F *hcorrmult56 = new TH2F("hcorrmult56","hcorrmult56",500,0,0.25,500,0,0.25);


   TH1F *hres45 = new TH1F("hres45","hres45",300,-5,1);
   TH1F *hres56 = new TH1F("hres56","hres56",300,-5,1);
   TH1F *hressum = new TH1F("hressum","hressum",300,-5,1);

   TH1F *hbin145 = new TH1F("hbin145","hbin145",300,-5,1);
   TH1F *hbin156 = new TH1F("hbin156","hbin156",300,-5,1);
   TH1F *hbin1sum = new TH1F("hbin1sum","hbin1sum",300,-5,1);
   TH1F *hbin245 = new TH1F("hbin245","hbin245",300,-5,1);
   TH1F *hbin256 = new TH1F("hbin256","hbin256",300,-5,1);
   TH1F *hbin2sum = new TH1F("hbin2sum","hbin2sum",300,-5,1);
   TH1F *hbin345 = new TH1F("hbin345","hbin345",300,-5,1);
   TH1F *hbin356 = new TH1F("hbin356","hbin356",300,-5,1);
   TH1F *hbin3sum = new TH1F("hbin3sum","hbin3sum",300,-5,1);
   TH1F *hbin445 = new TH1F("hbin445","hbin445",300,-5,1);
   TH1F *hbin456 = new TH1F("hbin456","hbin456",300,-5,1);
   TH1F *hbin4sum = new TH1F("hbin4sum","hbin4sum",300,-5,1);

   TH1F *hbin3mult45 = new TH1F("hbin3mult45","hbin3mult45",300,-5,1);
   TH1F *hbin3mult56 = new TH1F("hbin3mult56","hbin3mult56",300,-5,1);
   TH1F *hbin4mult45 = new TH1F("hbin4mult45","hbin4mult45",300,-5,1);
   TH1F *hbin4mult56 = new TH1F("hbin4mult56","hbin4mult56",300,-5,1);


   TH1F *hres45mult = new TH1F("hres45mult","hres45mult",300,-5,1);
   TH1F *hres56mult = new TH1F("hres56mult","hres56mult",300,-5,1);
   TH1F *hressummult = new TH1F("hressummult","hressummult",300,-5,1);

   TH1F *hn45220 = new TH1F("hn45220","hn45220",10,0,10);
   TH1F *hn56220 = new TH1F("hn56220","hn56220",10,0,10);
   TH1F *hn45mult = new TH1F("hn45mult","hn45mult",10,0,10);
   TH1F *hn56mult = new TH1F("hn56mult","hn56mult",10,0,10);

   TH1F *hxi45mult = new TH1F("hxi45mult","hxi45mult",100,0,0.25);
   TH1F *hxi56mult = new TH1F("hxi56mult","hxi56mult",100,0,0.25);
   TH1F *hxangle45mult = new TH1F("hxangle45mult","hxangle45mult",50,120,170);
   TH1F *hxangle56mult = new TH1F("hxangle56mult","hxangle56mult",50,120,170);
   TH1F *hxi45single = new TH1F("hxi45single","hxi45single",100,0,0.25);
   TH1F *hxi56single = new TH1F("hxi56single","hxi56single",100,0,0.25);
   TH1F *hxangle45single = new TH1F("hxangle45single","hxangle45single",50,120,170);
   TH1F *hxangle56single = new TH1F("hxangle56single","hxangle56single",50,120,170);

   TH1F *hmresmulti = new TH1F("hmresmulti","hmresmulti",500,-15,5);
   TH1F *hmressingle = new TH1F("hmressingle","hmressingle",500,-15,5);
   TH1F *hmresmixed = new TH1F("hmresmixed","hmresmixed",500,-15,5);
   TH2F *hmcorrmulti = new TH2F("hmcorrmulti","hmcorrmulti",250,0,2500,250,0,2500);
   TH2F *hmcorrmixed = new TH2F("hmcorrmixed","hmcorrmixed",250,0,2500,250,0,2500);
   TH2F *hmcorrsingle = new TH2F("hmcorrsingle","hmcorrsingle",250,0,2500,250,0,2500);
   TH1F *hmmumu = new TH1F("hmmumu","hmmumu",250,0,2500);
   TH2F *hycorrmulti = new TH2F("hycorrmulti","hycorrmulti",250,-2.5,2.5,250,-2.5,2.5);
   TH2F *hycorrsingle = new TH2F("hycorrsingle","hycorrsingle",250,-2.5,2.5,250,-2.5,2.5);
   TH2F *hycorrmixed = new TH2F("hycorrmixed","hycorrmixed",250,-2.5,2.5,250,-2.5,2.5);
   TH2F *hdmdysingle = new TH2F("hdmdysingle","hdmdysingle",1000,-500,500,200,-2,2);

   TH1F *hmresycutsingle = new TH1F("hmresycutsingle","hmresycutsingle",500,-15,5);
   TH1F *hmresycutmulti = new TH1F("hmresycutmulti","hmresycutmulti",500,-15,5);
   TH1F *hmresycutmixed = new TH1F("hmresycutmixed","hmresycutmixed",500,-15,5);

   TH1F *hpzmumumultmatch45 = new TH1F("hpzmumumultmatch45","hpzmumumultmatch45",500,-2000,2000);
   TH1F *hpzmumusinglematch45 = new TH1F("hpzmumusinglematch45","hpzmumusinglematch45",500,-2000,2000);
   TH1F *hpzmumumultmatch56 = new TH1F("hpzmumumultmatch56","hpzmumumultmatch56",500,-2000,2000);
   TH1F *hpzmumusinglematch56 = new TH1F("hpzmumusinglematch56","hpzmumusinglematch56",500,-2000,2000);
   TH1F *hpzmumumultantimatch45 = new TH1F("hpzmumumultantimatch45","hpzmumumultantimatch45",500,-2000,2000);
   TH1F *hpzmumusingleantimatch45 = new TH1F("hpzmumusingleantimatch45","hpzmumusingleantimatch45",500,-2000,2000);
   TH1F *hpzmumumultantimatch56 = new TH1F("hpzmumumultantimatch56","hpzmumumultantimatch56",500,-2000,2000);
   TH1F *hpzmumusingleantimatch56 = new TH1F("hpzmumusingleantimatch56","hpzmumusingleantimatch56",500,-2000,2000);

   TH2F *hxivst45 = new TH2F("hxivst45","hxivst45",500,0,0.25,100,0,5);
   TH2F *hxivst56 = new TH2F("hxivst56","hxivst56",500,0,0.25,100,0,5);

   ofstream ofs("TextOutputSingleMultiCorr.txt");

   ofs << "Run,LS,Event,Arm,xing angle,xi(p),xi(mumu),ximatch,t,theta*x,theta*y" << std::endl;

   Int_t nent = fChain->GetEntries();

   TLorentzVector mu1,mu2,mumu;

   Int_t usemultitracks = 0;
   Int_t ismc = 1;
   Int_t issideband = 0;
   Int_t year = 2018;
   Int_t usenear = 0;

   /*
    * Basic options for the type of selction to apply
    */
   usemultitracks = multi; // Allow only 1 proton per arm, or multiple protons
   ismc = mc; // Run on data or MC
   issideband = sb; // Select signal or "sideband" in acoplanarity+multiplicity
   year = yr; // Year
   usenear = nearfar; // For SingleRP, choose near(210) or far(220)

   Int_t id45 = 23;
   Int_t id56 = 123;
   if(usenear == 1)
     {
       id45 = 3;
       id56 = 103;
     }

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      if(jentry % 10000 == 0)
	std::cout << "Entry " << jentry << "/" << nent << std::endl;

      Float_t theacop = 1-fabs(Pair_dphi[0])/3.14159;

      // Basic selection cuts

      if(Pair_mass[0]<110) // M(mumu) > 110GeV
	continue;
      if(MuonCand_pt[0]<50 || MuonCand_pt[1]<50) // pT(mu) > 50 GeV
        continue;
      if(MuonCand_istight[0]!=1 || MuonCand_istight[1]!=1) // Tight muon ID
        continue;
      if(MuonCand_charge[0] == MuonCand_charge[1]) // Opposite charge muons
        continue;
      if((theacop > 0.009) && (issideband == 0)) // Dimuon acoplanarity <=0.009 for signal
	continue;
      if((theacop <= 0.009 || theacop > 0.1) && (issideband == 1)) // Acoplanarity 0.009-0.1 for sideband
        continue;
      if((Pair_extratracks0p5mm[0] > 0) && (issideband == 0)) // 0 extra tracks withing 0.5mm of dimuon vertex for signal
	continue;
      if((Pair_extratracks0p5mm[0] < 5 || Pair_extratracks0p5mm[0] > 10) && (issideband == 1)) // 5-10 extra tracks for sideband
        continue;

      mu1.SetPtEtaPhiE(MuonCand_pt[0],MuonCand_eta[0],MuonCand_phi[0],MuonCand_e[0]);
      mu2.SetPtEtaPhiE(MuonCand_pt[1],MuonCand_eta[1],MuonCand_phi[1],MuonCand_e[1]);
      mumu = mu1+mu2;


      hacop->Fill(1-fabs(Pair_dphi[0])/3.14159);

      Float_t mumuxisol1 = (1.0/13000.0)*((MuonCand_pt[0]*TMath::Exp(MuonCand_eta[0]) + MuonCand_pt[1]*TMath::Exp(MuonCand_eta[1])));
      Float_t mumuxisol2 = (1.0/13000.0)*((MuonCand_pt[0]*TMath::Exp(-1*MuonCand_eta[0]) + MuonCand_pt[1]*TMath::Exp(-1*MuonCand_eta[1])));

      Float_t protxi45 = 0.0;
      Float_t protxi56 = 0.0;
      Float_t protxi45single = 0.0;
      Float_t protxi56single = 0.0;
      Float_t protxi45multi = 0.0;
      Float_t protxi56multi = 0.0;
      Float_t protx45multi210 = 0.0;
      Float_t proty45multi210 = 0.0;
      Float_t protx45multi220 = 0.0;
      Float_t proty45multi220 = 0.0;
      Float_t protx56multi210 = 0.0;
      Float_t proty56multi210 = 0.0;
      Float_t protx56multi220 = 0.0;
      Float_t proty56multi220 = 0.0;
      Float_t mumuxi45 = 0.0;
      Float_t mumuxi56 = 0.0;
      Float_t prott45 = 0.0;
      Float_t prott56 = 0.0;
      Float_t protthx45 = 0.0;
      Float_t protthy45 = 0.0;
      Float_t protthx56 = 0.0;
      Float_t protthy56 = 0.0;
      Int_t ntrk45 = 0;
      Int_t ntrk56 = 0;

      Int_t ncountpixel45 = 0;
      Int_t ncountpixel56 = 0;
      Int_t ncountmulti45 = 0;
      Int_t ncountmulti56 = 0;

      // Count the number of single RP protons per pot
      for(Int_t p = 0; p < nRecoProtCand; p++)
        {
          if(ProtCand_rpid1[p] == id45 && ProtCand_ismultirp[p]==0)
            ncountpixel45++;
          if(ProtCand_rpid1[p] == id56 && ProtCand_ismultirp[p]==0)
            ncountpixel56++;
          if(ProtCand_ismultirp[p]==1 && ProtCand_arm[p]==0)
            ncountmulti45++;
          if(ProtCand_ismultirp[p]==1 && ProtCand_arm[p]==1)
            ncountmulti56++;
        }

      hn45220->Fill(ncountpixel45);
      hn56220->Fill(ncountpixel56);
      hn45mult->Fill(ncountmulti45);
      hn56mult->Fill(ncountmulti56);

      // Now do the actual analysis loop
      for(Int_t p = 0; p < nRecoProtCand; p++)
        {
	  // Sector 45, singleRP
          if(ProtCand_rpid1[p] == id45 && ProtCand_ismultirp[p]==0 && ntrk45 < 1 && (ncountpixel45==1 || usemultitracks==1))
            {
              protxi45 = ProtCand_xi[p];
              protxi45single = protxi45;
              hcorr45->Fill(protxi45,mumuxisol1);
              hres45->Fill(1 - (protxi45/mumuxisol1));
              hressum->Fill(1 - (protxi45/mumuxisol1));

              if(mumuxisol1 >= 0.02 && mumuxisol1 < 0.03)
                hbin145->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.03 && mumuxisol1 < 0.04)
                hbin245->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.04 && mumuxisol1 < 0.06)
                hbin345->Fill(1 - (protxi45/mumuxisol1));
              if(mumuxisol1 >= 0.06)
                hbin445->Fill(1 - (protxi45/mumuxisol1));

              if(fabs(1 - (protxi45/mumuxisol1)) < 0.25)
                {
                  hxi45single->Fill(protxi45);
                  hpzmumusinglematch45->Fill(mumu.Pz());
                }
              else
                hpzmumusingleantimatch45->Fill(mumu.Pz());

	      std::cout << "Run, Lumi, Event = " << Run << ", " << LumiSection << ", " << EventNum << std::endl;
	      std::cout << "\t xi(RP=3) = " << protxi45 << std::endl;

              //              ofs << Run << " " << LumiSection << " " << EventNum << " Single 23 " << protxi45 << " " << (1 - (protxi45/mumuxisol1)) << std::endl;              
              ntrk45++;
            }
          // Sector 56, singleRP                                                                                                                                                 
          if(ProtCand_rpid1[p] == id56 && ProtCand_ismultirp[p]==0 && ntrk56 < 1 && (ncountpixel56==1 || usemultitracks==1))
            {
              protxi56 = ProtCand_xi[p];
              protxi56single = protxi56;
              hcorr56->Fill(protxi56,mumuxisol2);
              hres56->Fill(1 - (protxi56/mumuxisol2));
              hressum->Fill(1 - (protxi56/mumuxisol2));

              if(mumuxisol2 >= 0.02 && mumuxisol2 < 0.03)
                hbin156->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.03 && mumuxisol2 < 0.04)
                hbin256->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.04 && mumuxisol2 < 0.06)
                hbin356->Fill(1 - (protxi56/mumuxisol2));
              if(mumuxisol2 >= 0.06)
                hbin456->Fill(1 - (protxi56/mumuxisol2));

              if(fabs(1 - (protxi56/mumuxisol2)) < 0.25)
                {
                  hxi56single->Fill(protxi56);
                  hpzmumusinglematch56->Fill(mumu.Pz());
                }
              else
                hpzmumusingleantimatch56->Fill(mumu.Pz());

	      std::cout << "Run, Lumi, Event = " << Run << " " << LumiSection << " " << EventNum << std::endl;
	      std::cout << "\t xi(RP=103) = " <<protxi56 << std::endl;

              //              ofs << Run << " " << LumiSection << " " << EventNum << " Single 123 " << protxi56 << " " << (1 - (protxi56/mumuxisol2)) << std::endl;             
              ntrk56++;
            }
          // Sector 45, multiRP                                                                                                                                                   
          if(ProtCand_arm[p] == 0 && ProtCand_ismultirp[p] == 1 && (ncountmulti45==1 || usemultitracks==1))
            {
              protxi45 = ProtCand_xi[p];
              prott45 = -1.0*ProtCand_t[p];
              protthx45 = ProtCand_ThX[p];
              protthy45 = ProtCand_ThY[p];
              protxi45multi = protxi45;

              // Tracks for eff. correction                                                                                                                                      
              protx45multi220 = ProtCand_trackx1[p];
              proty45multi220 = ProtCand_tracky1[p];
              protx45multi210 = ProtCand_trackx2[p];
              proty45multi210 = ProtCand_tracky2[p];

	      hres45mult->Fill(1 - (protxi45/mumuxisol1));
	      hressummult->Fill(1 - (protxi45/mumuxisol1));
	      hcorrmult45->Fill(protxi45,mumuxisol1);

	      if(mumuxisol1 >= 0.04 && mumuxisol1 < 0.06)
		hbin3mult45->Fill(1 - (protxi45/mumuxisol1));
	      if(mumuxisol1 >= 0.06)
		hbin4mult45->Fill(1 - (protxi45/mumuxisol1));

	      std::cout << "Run, Lumi, Event = " << Run << ", " << LumiSection << ", " << EventNum << std::endl;
	      std::cout << "\t xi(multi arm=45) = " << protxi45 << std::endl;

	      if(fabs(1 - (protxi45/mumuxisol1)) < 0.10 && (mumuxisol1 >= 0.04))
		{
		  hxi45mult->Fill(protxi45);
		  hxangle45mult->Fill(CrossingAngle);
		  hpzmumumultmatch45->Fill(mumu.Pz());
		  hxivst45->Fill(protxi45,prott45);

		  ofs << Run << "," << LumiSection << "," << EventNum << ",45," << CrossingAngle << "," << protxi45 << "," << mumuxisol1 << ","
		      << 1 - (protxi45/mumuxisol1) << "," << prott45
		      << ", " << protthx45 << ", " << protthy45 << ", " << ", " 
		      << std::endl;
		}
	      else
		hpzmumumultantimatch45->Fill(mumu.Pz());
	    }
	  // Sector 56, multiRP
	  if(ProtCand_arm[p] == 1 && ProtCand_ismultirp[p] == 1 && (ncountmulti56==1 || usemultitracks==1))
	    {
	      protxi56 = ProtCand_xi[p];
	      prott56 = -1.0*ProtCand_t[p];
	      protthx56 = ProtCand_ThX[p];
	      protthy56 = ProtCand_ThY[p];
	      protxi56multi = protxi56;
	      
	      // Tracks for eff. correction                                                                                                                                      
	      protx56multi220 = ProtCand_trackx1[p];
	      proty56multi220 = ProtCand_tracky1[p];
	      protx56multi210 = ProtCand_trackx2[p];
	      proty56multi210 = ProtCand_tracky2[p];
	      
	      hres56mult->Fill(1 - (protxi56/mumuxisol2));
	      hressummult->Fill(1 - (protxi56/mumuxisol2));
	      hcorrmult56->Fill(protxi56,mumuxisol2);
	      
	      if(mumuxisol2 >= 0.04 && mumuxisol2 < 0.06)
		hbin3mult56->Fill(1 - (protxi56/mumuxisol2));
	      if(mumuxisol2 >= 0.06)
		hbin4mult56->Fill(1 - (protxi56/mumuxisol2));
	      
	      
	      std::cout << "Run, Lumi, Event = " << Run << ", " << LumiSection << ", " << EventNum << std::endl;
	      std::cout << "\t xi(multi arm=56) = " << protxi56 << std::endl;
	      
	      if((fabs(1 - (protxi56/mumuxisol2)) > -0.05) && (fabs(1 - (protxi56/mumuxisol2)) < 0.20) && (mumuxisol2 >= 0.04))
		{
		  hxi56mult->Fill(protxi56);
		  hxangle56mult->Fill(CrossingAngle);
		  hpzmumumultmatch56->Fill(mumu.Pz());
		  hxivst56->Fill(protxi56,prott56);
		  
		  ofs << Run << "," << LumiSection << "," << EventNum << ",56," << CrossingAngle << "," << protxi56 << "," << mumuxisol2 << ","
		      << 1 - (protxi56/mumuxisol2) << "," << prott56
		      << ", " << protthx56 << ", " << protthy56 << ", " << ", " 
		      << std::endl;
		  
		}
	      else
		hpzmumumultantimatch56->Fill(mumu.Pz());
	    }
	}
   }

   // Textfile output
   ofs.close();
          
   // ROOT histogram output
   TFile *fx;
   fx = new TFile("DimuonHistograms.root","RECREATE");

   // Write all histograms
   hacop->Write();
   hcorr45->Write();
   hres45->Write();
   hcorr56->Write();
   hres56->Write();
   hressum->Write();
   hres45mult->Write();
   hres56mult->Write();
   hressummult->Write();
   
   hbin145->Write();
   hbin245->Write();
   hbin345->Write();
   hbin445->Write();

   hbin156->Write();
   hbin256->Write();
   hbin356->Write();
   hbin456->Write();

   hbin3mult45->Write();
   hbin3mult56->Write();
   hbin4mult45->Write();
   hbin4mult56->Write();

   hn45220->Write();
   hn56220->Write();
   hn45mult->Write();
   hn56mult->Write();

   hcorrmult45->Write();
   hcorrmult56->Write();

   hxi45mult->Write();
   hxi56mult->Write();
   hxangle45mult->Write();
   hxangle56mult->Write();
   hxi45single->Write();
   hxi56single->Write();

   hmresmulti->Write();
   hmresmixed->Write();
   hmressingle->Write();
   hmcorrsingle->Write();
   hmcorrmulti->Write();
   hmcorrmixed->Write();

   hmmumu->Write();
   hycorrsingle->Write();
   hycorrmulti->Write();
   hycorrmixed->Write();
   hdmdysingle->Write();
   hmresycutmulti->Write();
   hmresycutmixed->Write();
   hmresycutsingle->Write();
   
   hpzmumusinglematch45->Write();
   hpzmumusinglematch56->Write();
   hpzmumumultmatch45->Write();
   hpzmumumultmatch56->Write();
   hpzmumusingleantimatch45->Write();
   hpzmumusingleantimatch56->Write();
   hpzmumumultantimatch45->Write();
   hpzmumumultantimatch56->Write();

   hxivst45->Write();
   hxivst56->Write();

   fx->Write();

}
