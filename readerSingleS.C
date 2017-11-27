//  ITU CMS Group Analysis Framework
//undfined DelphMT2W isoTracksP DelphMT
//         MT2W0                MT2W0
#include "NtupleTools3.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TFile.h>
#include <TBranch.h>

using namespace std;

void readerSingleS(TString list, TString outname,bool useW=true){

	bool D=false;// DelphMet
	
	TObjArray* arr = list.Tokenize(" ");
	int size=arr->GetEntries();
	if(size%2!=0) {
	  cout<<"unbalance file/weight list: "<<list<<endl;
	  exit(0);
	}
	vector<TString> files;
	vector<double> weights;
	for(int i=0;i<size;i+=2){
	  files.push_back( arr->At(i)->GetName() );
	  weights.push_back( atof( arr->At(i+1)->GetName() ) );
	}
	
	//---------Define tree variables-------------------------------
	vector<double> HT_CUT;
	vector<double> MT_CUT;
	vector<double> MET_CUT;
	vector<double> MT2W_CUT;
		
	//--------------------------------------------------------------
	TTree *TreeB= new TTree("TreeB","Background Tree");
	TTree *TreeS= new TTree("TreeS","Signal Tree");	

	TreeB->Branch("HT",&HT_CUT);
	TreeB->Branch("MET",&MET_CUT);
	TreeB->Branch("MT_CUT",& MT_CUT);
	TreeB->Branch("MT2W",&MT2W_CUT);
	
	TreeS->Branch("HT",&HT_CUT);
	TreeS->Branch("MET",&MET_CUT);
	TreeS->Branch("MT_CUT",&MT_CUT);
	TreeS->Branch("MT2W",&MT2W_CUT);
	

        //at 1lep 4jets 1b
        TH1F* aDphi      = new TH1F("aDphi","#Delta#phi",100,0,4);
        TH1F* aDelphDphi = new TH1F("aDelphDphi","Delphes #Delta#phi",100,0,4);

        TH1F* aMET       = new TH1F("aMET","MET",100,0,3000);
        TH1F* aDelphMET  = new TH1F("aDelphMET","Delphes MET",100,0,3000);

        TH1F* aMT        = new TH1F("aMT","MT",100,0,500);
        TH1F* aDelphMT   = new TH1F("aDelphMT","Delphes MT",100,0,500);
        TH1F* aMT2       = new TH1F("aMT2","MT",100,0,500);
        TH1F* aDelphMT2  = new TH1F("aDelphMT2","Delphes MT",100,0,500);


        TH1F* aMT2W       = new TH1F("aMT2W","MT2W 1,2b only",100,0,500);
        TH1F* aDelphMT2W  = new TH1F("aDelphMT2W","Delphes MT2W 1,2b only",100,0,500);
        TH1F* aMT2W2      = new TH1F("aMT2W2","MT2W 1,2b only",100,0,500);
        TH1F* aDelphMT2W2 = new TH1F("aDelphMT2W2","Delphes MT2W 1,2b only",100,0,500);

        TH1F* aHT        = new TH1F("aHT","HT40",100,0,4000);
        TH1F* aLepPt     = new TH1F("aLepPt","single lep Pt",100,0,1000);
        TH1F* aLepEta     = new TH1F("aLepEta","single lep Eta",100,-4.5,4.5);
        TH1F* a1JetPt    = new TH1F("a1JetPt","highest jet Pt",100,0,2000);
        TH1F* a1bJetPt   = new TH1F("a1bJetPt","highest bjet Pt",100,0,2000);
        TH1F* an1bJetPt   = new TH1F("an1bJetPt","2bd and more highest bjets Pt",100,0,2000);

	TH1F* hMET   = new TH1F("hMET","MET",100,0,3000);
	TH1F* hMT    = new TH1F("hMT","MT",100,0,300);
	TH1F* hMT2Wpre  = new TH1F("hMT2Wpre","MT2W w/o MET,MT req.",100,0,500);
	TH1F* hMT2W  = new TH1F("hMT2W","MT2W",100,0,500);
	TH1F* hMETMeff  = new TH1F("hMETMeff","MET/(MET+HT)",100,0,1);
	TH1F* hHT    = new TH1F("hHT","HT40",100,0,5000);
	TH1F* hDphi  = new TH1F("hDphi","#Delta#phi",100,0,4);
	TH1F* hAllHT = new TH1F("hAllHT","before cuts HT40",100,0,5000);
	TH1F* hAllLepPt = new TH1F("hAllLepPt","singlw lep Pt - all",100,0,2000);
	TH1F* hLepPtJb = new TH1F("hLepPtJ","single lep Pt - jet req",100,0,2000);
	TH1F* hLepPtM = new TH1F("hLepPtM","single lep Pt - jet + MET",100,0,2000);
	TH1F* hLepPtMM = new TH1F("hLepPtMM","single lep Pt - jet + MET +Mt+MT2W",100,0,2000);
	TH1F* hLepEtaJb = new TH1F("hLepEtaJ","single lep eta - jet req",100,-4.5,4.5);
	TH1F* hLepEtaM = new TH1F("hLepEtaM","single lep eta - jet + MET",100,-4.5,4.5);
	TH1F* hLepEtaMM = new TH1F("hLepEtaMM","single lep eta - jet + MET +Mt+MT2W",100,-4.5,4.5);


        TH1F* abasym   = new TH1F("abasym","(bjetsPt-lepPt)/(bjetsPt+lepPt)",100,-1,1);
        TH1F* hbasym   = new TH1F("hbasym","(bjetsPt-lepPt)/(bjetsPt+lepPt)",100,-1,1);
        TH1F* aC      = new TH1F("aC","centrality",100,0,5);
        TH1F* hC      = new TH1F("hC","centrality",100,0,5);
	TH1F* atop = new TH1F("atop","topness",100,-21,21);
	TH1F* atop0 = new TH1F("atop0","topness Delph MET",100,-21,21);
	TH1F* htop = new TH1F("htop","topness",100,-21,21);
	TH1F* htop0 = new TH1F("htop0","topness Delph MET",100,-21,21);
	TH1F* htHT    = new TH1F("htHT","topness HT40",100,0,5000);
	TH1F* htDphi  = new TH1F("htDphi","topness #Delta#phi",100,0,4);

	TH1I* hNel = new TH1I("hNel","N good(30) el",10,0,10);
	TH1I* hNmu = new TH1I("hNmu","N good(30) mu",10,0,10);
	TH1I* hNlep = new TH1I("hNlep","N good(30) mu+el",10,0,10);
	TH1I* hNtj = new TH1I("hNtj","N good(40) jets",15,0,15);
	TH1I* hNtj_cut1 = new TH1I("hNtj_cut1","N good(40) jets, Njet > 4",15,0,15);
        TH1I* hNbjet = new TH1I("hNbjet","N b jets",15,0,15);

	EasyChain* tree = new EasyChain("delphTree");
	
	for(unsigned i=0;i<files.size();i++){
       		tree->AddSmartW(files[i],weights[i]);
		cout<<"add: "<<files[i]<<" "<<weights[i]<<endl;
	}

	int Nevents=tree->GetEntries();
	cout<<">>>>>>>>>>>>>>>>>>>>>>> total number of events:\t" << Nevents <<endl;

	// CutFlow variables
	const int CutNumb = 9;
	const char * CutList[CutNumb] = {"#Events","sngl. lep.",
					 "nJets >= 4","bjets >= 2",
                                         "MET>500","dphi12 > 0.5",
                                         "MT>120",
                                         "--MT2W>150",
                                         "HT>500"
                                         };
	
	double CFCounter[CutNumb];
	int   iCFCounter[CutNumb];
	
	for (int i=0;i < CutNumb; i++){
	  CFCounter[i] = 0;
	  iCFCounter[i] = 0;
	}
	TH1D* CutFlow= new TH1D("CutFlow","Cut Flow",CutNumb,0.5,CutNumb+0.5);
	// label bins
        for(int cj = 0; cj < CutNumb; cj++)
	  CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj]);
	HT_CUT.clear();
        MET_CUT.clear();
	MT_CUT.clear();
        MT2W_CUT.clear();           
	

	for(int entry=0; entry < Nevents; entry+=1){
//		progress();

		progressT();
		double fw = tree->GetEntryW(entry); // the applied with AddSmartW for the current file/dir

		double EvWeight = 1;
		if(useW) EvWeight = tree->Get(EvWeight,"EventWeight");
		EvWeight *= fw * 1000;

//		double HT   = tree->Get(HT,"HT"); 
		double HT   = tree->Get(HT,"HT"); 

		// 0. CF presel
	        CFCounter[0]+= EvWeight;
		iCFCounter[0]++;
		
		std::vector<TLorentzVector> &Electrons = tree->Get(&Electrons,"Electrons");
		std::vector<TLorentzVector> &Muons = tree->Get(&Muons,"Muons");
		std::vector<int> &ElectronCh = tree->Get(&ElectronCh,"ElectronCh");
		std::vector<int> &MuonCh     = tree->Get(&MuonCh,"MuonCh");

		// 1. Lepton veto
		// etamin in ntupler 2.5, ptmin 10
		int    Nel_tight = 0;
		int    Nmu_tight = 0;
		int    Nel_loose = Electrons.size();
		int    Nmu_loose = Muons.size();
		// 
		double lepPt=0;
		double lepEta=0;
		TLorentzVector* Lep=0;
		for(unsigned i=0;i<Nmu_loose;++i) 
			if(fabs(Muons[i].Eta())<=2.4&&Muons[i].Pt()>30 ){
				Nmu_tight++;
				lepPt=Muons[0].Pt();
				lepEta=Muons[0].Eta();
				Lep=&Muons[0];
                         	}
		for(unsigned i=0;i<Nel_loose;++i) 
			if(fabs(Electrons[i].Eta())<=2.4&&Electrons[i].Pt()>30 ){
				Nel_tight++;
				lepPt=Electrons[0].Pt();
				lepEta=Electrons[0].Eta();
				Lep=&Electrons[0];
			}
		
		// exactly 1 hard lepton, no other loose
		if(Nel_tight+Nmu_tight != 1) continue;
                

		// 2. nJets >= 3,4
		vector<TLorentzVector> &Jets = tree->Get(&Jets,"Jets");
		int Njet_loose = tree->Get(Njet_loose,"Njet");
		int Njet_tight = 0;
		for(int i = 0;i  < Njet_loose; i++)
			if(Jets[i].Pt() > 40) Njet_tight++;
                          
	        CFCounter[1]+= EvWeight;
		iCFCounter[1]++;



		// distributions at good object selections ----------------------------------------------
                
                hNmu->Fill(Nmu_tight,EvWeight);
                hNmu->Fill(Nel_tight,EvWeight);
                hNlep->Fill(Nel_tight+Nmu_tight,EvWeight);
                hAllLepPt->Fill(lepPt,EvWeight);
                hAllHT->Fill(HT,EvWeight);
                hNtj->Fill(Njet_tight,EvWeight);

               //----------------------------------------------------------------------------------------


		if(Njet_tight < 4) continue;
		 CFCounter[2]+= EvWeight;
		iCFCounter[2]++;
	
                hNtj_cut1->Fill(Njet_tight,EvWeight);
	

		// 3. Btag cut
		double bjetpt=0;
		int Nbjet = 0;
		bool hardB = false;
		vector<int> &JetB = tree->Get(&JetB,"JetB");		
		for(unsigned i=0;i<Jets.size();i++) {
			if(JetB[i]>0&&Jets[i].Pt() > 40) {
				if(Nbjet==0) {
					a1bJetPt->Fill(Jets[i].Pt());
					bjetpt=Jets[i].Pt();
				}
				else         an1bJetPt->Fill(Jets[i].Pt());
				Nbjet++;
			}
			if(JetB[i]>0&&Jets[i].Pt() > 250) hardB=true;
		}
		//if(Nbjet<1||Nbjet>2) continue;		
		if(Nbjet<1) continue; 
                hNbjet->Fill(Nbjet,EvWeight);
                hLepPtJb->Fill(lepPt,EvWeight);
                hLepEtaJb->Fill(lepEta,EvWeight);
         
                CFCounter[3]+= EvWeight;
		iCFCounter[3]++;
		
                // distributions at 1lep4jet1b
                
                double basym = (bjetpt-lepPt)/(bjetpt+lepPt);
                abasym->Fill(basym);
                
                vector<double> &JetMETdPhi = tree->Get(&JetMETdPhi,"JetMETdPhi");
                double dPhi = TMath::Min(JetMETdPhi[0],JetMETdPhi[1]);
                aDphi->Fill(dPhi,EvWeight);

                double DelphMET_Phi  = tree->Get(DelphMET_Phi,"DelphMET_Phi");
                double DelphdPhi = TMath::Min(acos(cos(Jets[0].Phi()-DelphMET_Phi)),acos(cos(Jets[1].Phi()-DelphMET_Phi)));
                aDelphDphi->Fill(DelphdPhi,EvWeight);
		
                double MET  = tree->Get(MET,"MET");
                aMET->Fill(MET,EvWeight);
                double DMET  = tree->Get(DMET,"DelphMET");
                aDelphMET->Fill(DMET,EvWeight);
                
                double MT2W  = tree->Get(MT2W,"MT2W");
                aMT2W2->Fill(MT2W,EvWeight);
                //if(MT2W>0&&MT2W<499)aMT2W->Fill(MT2W,EvWeight);
                double DMT2W  = tree->Get(DMT2W,"MT2W0");
                aDelphMT2W2->Fill(DMT2W,EvWeight);

                double MT  = tree->Get(MT,"MT");
	        aMT2->Fill(MT,EvWeight);
                if(MT!=0)aMT->Fill(MT,EvWeight);
                double DMT  = tree->Get(DMT,"MT0");
                aDelphMT2->Fill(DMT,EvWeight);
                if(DMT!=0)aDelphMT->Fill(DMT,EvWeight);

                //aHT->Fill(HT,EvWeight);
                aLepPt->Fill(lepPt,EvWeight);
                aLepEta->Fill(lepEta,EvWeight);
                if(Jets[0].Pt()>40.) a1JetPt->Fill(Jets[0].Pt(),EvWeight);

                //double top   = tree->Get(top,"top");
                //double top0  = tree->Get(top0,"top0");

		//if(top!=0)atop->Fill(top,EvWeight);
		//if(top0!=0)atop0->Fill(top0,EvWeight);

		// Centrality
                double MET_Phi  = tree->Get(MET_Phi,"MET_Phi");
		double METx,METy;
		if(D){
		METx = cos(DelphMET_Phi)*DMET;
	        METy = sin(DelphMET_Phi)*DMET;
		}else{
		METx = cos(MET_Phi)*MET;
	        METy = sin(MET_Phi)*MET;
		}
                TLorentzVector a=Jets[0]+(*Lep);
                TLorentzVector b=Jets[1]+(*Lep);
                TLorentzVector c=b;
                b.SetPx(b.Px()+METx);
                b.SetPy(b.Py()+METy);
                double deta1=fabs(a.Eta()-b.Eta());
                a.SetPx(b.Px()+METx);
                a.SetPy(b.Py()+METy);
                double deta2=fabs(a.Eta()-c.Eta());
                double C=TMath::Max(deta1,deta2);
		aC->Fill(C,EvWeight);


		//  MET/Meff
		double Meff = MET + HT;
		if(D)  Meff = DMET + HT;
		if(D) { if(DMET>0)hMETMeff->Fill(DMET/Meff,EvWeight); } // don't remove the curly bracket!
		else  if(MET>0)hMETMeff->Fill(MET/Meff,EvWeight);
		
		// 4. MET cut
		if(D){
		 hMET->Fill(DMET,EvWeight);
		 if(DMET<250) continue;
		}else{
		 hMET->Fill(MET,EvWeight);
		 if(MET<250) continue;
		}
		CFCounter[4]+= EvWeight;
		iCFCounter[4]++;
		

		// 6. dPhi  - had been inverted in first version
		if(D){
		 hDphi->Fill(DelphdPhi,EvWeight);
		 if(DelphdPhi < 0.5) continue;
		}else{
		 hDphi->Fill(dPhi,EvWeight);
		 if(dPhi < 0.5) continue;
		}

		 CFCounter[5]+= EvWeight;
		iCFCounter[5]++;
		

		hLepPtM->Fill(lepPt,EvWeight);
		hLepEtaM->Fill(lepEta,EvWeight);

		// 5. MT cut
//		double MT  = tree->Get(MT,"MT");
		if(D){
		 hMT->Fill(DMT,EvWeight);
		 if(DMT<120) continue;
		}else{
		 hMT->Fill(MT,EvWeight);
		 if(MT<120) continue;
		}
		 CFCounter[6]+= EvWeight;
		iCFCounter[6]++;
                

                if(D){ if(DMT2W<170) continue; hMT2Wpre->Fill(DMT2W,EvWeight);}
                else if(MT2W<170) continue; hMT2Wpre->Fill(MT2W,EvWeight);
                
                CFCounter[7]+= EvWeight;
                iCFCounter[7]++;

                if(HT < 500) continue;
                aHT->Fill(HT,EvWeight);   

                CFCounter[8]+= EvWeight;
                iCFCounter[8]++;

		HT_CUT.push_back(HT);
		MET_CUT.push_back(MET);
		MT2W_CUT.push_back(MT2W); 
		MT_CUT.push_back(MT);

		



                                

		//TreeS->Fill();
	}
	// ^loop end^

        TreeS->Fill();

	ofstream tfile;
	if(D)tfile.open("SingleSDelphMET_"+outname+".txt");
	else tfile.open("SingleS_"+outname+".txt");
	tfile << "########################################" << endl;
	tfile << "Cut efficiency numbers:" << endl;
	for(int ci = 0; ci < CutNumb; ci++)
	{
		tfile << "After cut " << CutList[ci] <<  "\t\t\t"
		      << CFCounter[ci]  << "\t events left\t"<< iCFCounter[ci] <<" cnt\t"<< endl;
		CutFlow->SetBinContent(1+ci,CFCounter[ci]);
	}

	TFile * outf;
	if(D) outf  = new TFile("SingleSDelphMET_"+outname+"_his.root","RECREATE");
	else  outf  = new TFile("SingleS_"+outname+"_his.root","RECREATE");

	hDphi->Write();
	hHT->Write();
	hAllHT->Write();
	hMETMeff->Write();
	hMET->Write();
	hMT->Write();
	hMT2W->Write();
	hNel->Write();
	hNmu->Write();
	hNtj->Write();
        hNtj_cut1->Write();
	hNel->Write();
	hNmu->Write();
	hNlep->Write();
	hNbjet->Write();

	hAllLepPt->Write();
	hLepPtJb->Write();
	hLepPtM->Write();
	hLepPtMM->Write();
	hLepEtaJb->Write();
	hLepEtaM->Write();
	hLepEtaMM->Write();

	hMT2Wpre->Write();

        aDphi->Write();
        aDelphDphi->Write();

        aMET->Write();
        aDelphMET->Write();

        aMT->Write();
        aMT2->Write();
        aDelphMT->Write();
        aDelphMT2->Write();

        aMT2W->Write();
        aMT2W2->Write();
        aDelphMT2W->Write();
        aDelphMT2W2->Write();

        aHT->Write();
        aLepPt->Write();
        aLepEta->Write();
        a1JetPt->Write();
        a1bJetPt->Write();
        an1bJetPt->Write();

	//atop->Write();
	//atop0->Write();
	//htop->Write();
	//htop0->Write();
	htHT->Write();
	htDphi->Write();
	abasym->Write();
	hbasym->Write();
	aC->Write();
	hC->Write();    
	TreeS->Write();
	TreeB->Write();
//	for(int i=0;i<check.size();i++){
//		cout<<check[i]<<" "<<CutList[i]<<endl;
//	}
}
