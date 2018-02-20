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

void readerSingleS_cutflows(TString list, TString outname,bool useW=true){
    
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
    vector<double> VEC_AWH2;
    vector<double> VEC_AWH1;
    
    
    //--------------------------------------------------------------
    TTree *TreeB0 = new TTree("TreeB0","TTBar");
    TTree *TreeB1 = new TTree("TreeB1","Top Jets");
    TTree *TreeB2 = new TTree("TreeB2","Boson Jets");
    TTree *TreeB3 = new TTree("TreeB3","Diboson");
    TTree *TreeS  = new TTree("TreeS","Signal Tree");
    
    TreeB0->Branch("HT",&HT_CUT);
    TreeB0->Branch("MET",&MET_CUT);
    TreeB0->Branch("MT",&MT_CUT);
    TreeB0->Branch("MT2W",&MT2W_CUT);
    TreeB0->Branch("AWH2",&VEC_AWH2);
    TreeB0->Branch("AWH1",&VEC_AWH1);
    
    TreeB1->Branch("HT",&HT_CUT);
    TreeB1->Branch("MET",&MET_CUT);
    TreeB1->Branch("MT",&MT_CUT);
    TreeB1->Branch("MT2W",&MT2W_CUT);
    TreeB1->Branch("AWH2",&VEC_AWH2);
    TreeB1->Branch("AWH1",&VEC_AWH1);
    
    TreeB2->Branch("HT",&HT_CUT);
    TreeB2->Branch("MET",&MET_CUT);
    TreeB2->Branch("MT",&MT_CUT);
    TreeB2->Branch("MT2W",&MT2W_CUT);
    TreeB2->Branch("AWH2",&VEC_AWH2);
    TreeB2->Branch("AWH1",&VEC_AWH1);
    
    TreeB3->Branch("HT",&HT_CUT);
    TreeB3->Branch("MET",&MET_CUT);
    TreeB3->Branch("MT",&MT_CUT);
    TreeB3->Branch("MT2W",&MT2W_CUT);
    TreeB3->Branch("AWH2",&VEC_AWH2);
    TreeB3->Branch("AWH1",&VEC_AWH1);
    
    TreeS->Branch("HT",&HT_CUT);
    TreeS->Branch("MET",&MET_CUT);
    TreeS->Branch("MT",&MT_CUT);
    TreeS->Branch("MT2W",&MT2W_CUT);
    TreeS->Branch("AWH2",&VEC_AWH2);
    TreeS->Branch("AWH1",&VEC_AWH1);
    
    //Histos
    TH1F* hHT_0      = new TH1F("hHT_0","HT",100,0,4000);
    TH1F* hLepPt_0   = new TH1F("hLepPt_0","Single Lepton pT",100,0,2000);
    TH1F* hLepEta_0   = new TH1F("hLepEta_0","Single lep eta",100,-4.5,4.5);
    TH1F* hHT_1      = new TH1F("hHT_1","HT",100,0,4000);
    TH1F* hLepPt_1   = new TH1F("hLepPt_1","Single Lepton pT",100,0,2000);
    TH1F* hLepEta_1   = new TH1F("hLepEta_1","Single lep eta",100,-4.5,4.5);
    TH1F* hJetPt_1     = new TH1F("hJetPt_1","highest jet Pt",100,0,2000);
    TH1F* hMET_cut       = new TH1F("hMET_cut","MET",100,0,3000);
    TH1F* hDphi_cut      = new TH1F("hDphi_cut","#Delta#phi",100,0,3.2);
    TH1F* hLepPt_2   = new TH1F("hLepPt_2","Single Lepton pT",100,0,2000);
    TH1F* hLepEta_2   = new TH1F("hLepEta_2","Single lep eta",100,-4.5,4.5);
    TH1F* hMT_cut        = new TH1F("hMT_cut","MT",100,0,500);
    TH1F* hHT_cut      = new TH1F("hHT_cut","HT",100,0,4000);
    TH1F* hLepPt_3   = new TH1F("hLepPt_3","Single Lepton pT",100,0,2000);
    TH1F* hLepEta_3   = new TH1F("hLepEta_3","Single lep eta",100,-4.5,4.5);
    
    TH1I* hNel   = new TH1I("hNel","N good(30) el",10,0,10);
    TH1I* hNmu      = new TH1I("hNmu","N good(30) mu",10,0,10);
    TH1I* hNjet   = new TH1I("hNjet","N good(40) jets", 15, 0, 15);
    TH1I* hNjet_check   = new TH1I("hNjet_check","N good(40) jets, Njet > 4",15,0,15);
    TH1F* hmaxbJetPt    = new TH1F("hmaxbJetPt","Highest Bjet Pt",100,0,2000);
    TH1F* hmaxpbJetPt   = new TH1F("hmaxpbJetPt","2^{nd} Highest and More Bjets Pt",100,0,2000);
    TH1I* hNbjet   = new TH1I("hNbjet","N b jets",15,0,15);
    TH1F* hbasym     = new TH1F("hbasym","(bjetsPt-lepPt)/(bjetsPt+lepPt)",100,-1,1);
    TH1F* hDphi      = new TH1F("hDphi","#Delta#phi",100,0,3.2);
    TH1F* hDelphDphi = new TH1F("hDelphDphi","Delphes #Delta#phi",100,0,3.2);
    TH1F* hMET       = new TH1F("hMET","MET",100,0,3000);
    TH1F* hDelphMET  = new TH1F("hDelphMET","Delphes MET",100,0,3000);
    TH1F* hMT2W      = new TH1F("hMT2W","MT2W",100,0,600);
    TH1F* hDelphMT2W = new TH1F("hDelphMT2W","Delphes MT2W",100,0,500);
    TH1F* hMT        = new TH1F("hMT","MT",100,0,500);
    TH1F* hMT0        = new TH1F("hMT0","if MT != 0",100,0,500);
    TH1F* hDelphMT   = new TH1F("hDelphMT","Delphes MT",100,0,500);
    TH1F* hDelphMT0   = new TH1F("hDelphMT0"," if Delphes MT != 0",100,0,500);
    TH1F* hCentrality = new TH1F("hCentrality","centrality",100,0,5);
    TH1F* hMETMeff    = new TH1F("hMETMeff","MET/(MET+HT)",100,0,1);
    TH1F* hAW     = new TH1F("hAW","Meff/ #sqrt{HT}",100,0,100);
    TH1F* hAW1    = new TH1F("hAW1","MT2W/ #sqrt{HT}",100,0,40);
    TH1F* hAW2    = new TH1F("hAW2","MT2W/ #sqrt{HT+MET}",100,0,35);
    TH1F* hAW3    = new TH1F("hAW3","MT2W/ #sqrt{MT+MET}",10,0,60);
    TH1F* hAW4    = new TH1F("hAW4","MT2W/ #sqrt{MT+MET+HT}",10,0,60);
    TH1F* hMT2Wpre    = new TH1F("hMT2Wpre","MT2W w/o MET,MT req.",100,0,500);
    TH1F* hDenom      = new TH1F("hDenom","#sqrt{HT+MET}",100,0,60);
    TH1F* hsHT       = new TH1F("hsHT","#sqrt{HT}",100,0,50);
    TH1F* hsMET   = new TH1F("hsMET","#sqrt{MET}",100,0,1);
    //Histos
    
    EasyChain* tree      = new EasyChain("delphTree");
    
    for(unsigned i=0;i<files.size();i++){
        tree->AddSmartW(files[i],weights[i]);
        cout<<"add: "<<files[i]<<" "<<weights[i]<<endl;
    }
    
    int Nevents=tree->GetEntries();
    cout<<">>>>>>>>>>>>>>>>>>>>>>> total number of events:\t" << Nevents <<endl;
    
    // CutFlow variables
    const int CutNumb = 9; //number of cuts
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
    for(int cj = 0; cj < CutNumb; cj++)
        CutFlow->GetXaxis()->SetBinLabel(cj+1,CutList[cj]);
    
    HT_CUT.clear();
    MET_CUT.clear();
    MT_CUT.clear();
    MT2W_CUT.clear();
    VEC_AWH2.clear();
    VEC_AWH1.clear();
    
    for(int entry=0; entry < Nevents; entry+=1){
        progressT();
        double fw = tree->GetEntryW(entry); // the applied with AddSmartW for the current file/dir
        double EvWeight = 1.0;
        
        if(useW) EvWeight = tree->Get(EvWeight,"EventWeight");
        EvWeight *= fw * 1000 ;
        
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
        for(unsigned i=0; i<Nmu_loose;++i)
            if(fabs(Muons[i].Eta())<=2.4&&Muons[i].Pt()>30 ){
                Nmu_tight++;
                lepPt=Muons[0].Pt();
                lepEta=Muons[0].Eta();
                Lep=&Muons[0];
            }
        for(unsigned i=0; i<Nel_loose;++i)
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
            if(Jets[i].Pt() > 40){
                Njet_tight++;
            }
        
        CFCounter[1]+= EvWeight;
        iCFCounter[1]++;
        
        //Good Event Selection
        hHT_0->Fill(HT,EvWeight);
        hNel->Fill(Nel_tight,EvWeight);
        hNmu->Fill(Nmu_tight,EvWeight);
        hNjet->Fill(Njet_tight,EvWeight);
        hLepPt_0->Fill(lepPt,EvWeight);
        hLepEta_0->Fill(lepEta,EvWeight);
        //Good Event Selection
        
        if(Njet_tight < 4) continue;
        CFCounter[2]+= EvWeight;
        iCFCounter[2]++;
        
        hNjet_check->Fill(Njet_tight,EvWeight);
        
        // 3. Btag cut
        double bjetpt=0;
        int Nbjet = 0;
        bool hardB = false;
        vector<int> &JetB = tree->Get(&JetB,"JetB");
        for(unsigned i=0;i<Jets.size();i++) {
            if(JetB[i]>0&&Jets[i].Pt() > 40) {
                if(Nbjet==0) {
                    hmaxbJetPt->Fill(Jets[i].Pt());
                    bjetpt=Jets[i].Pt();
                }
                else hmaxpbJetPt->Fill(Jets[i].Pt());
                Nbjet++;
            }
            if(JetB[i]>0&&Jets[i].Pt() > 250) hardB=true;
        }
        
        //if(Nbjet<1||Nbjet>2) continue;
        if(Nbjet<1) continue;
        
        hNbjet->Fill(Nbjet,EvWeight);
        hLepPt_1->Fill(lepPt,EvWeight);
        hLepEta_1->Fill(lepEta,EvWeight);
        hHT_1->Fill(HT,EvWeight);
        if(Jets[0].Pt()>40.) hJetPt_1->Fill(Jets[0].Pt(),EvWeight);
        
        CFCounter[3]+= EvWeight;
        iCFCounter[3]++;
        
        double basym = (bjetpt-lepPt)/(bjetpt+lepPt);
        vector<double> &JetMETdPhi = tree->Get(&JetMETdPhi,"JetMETdPhi");
        double dPhi = TMath::Min(JetMETdPhi[0],JetMETdPhi[1]);
        double DelphMET_Phi  = tree->Get(DelphMET_Phi,"DelphMET_Phi");
        double DelphdPhi = TMath::Min(acos(cos(Jets[0].Phi()-DelphMET_Phi)),acos(cos(Jets[1].Phi()-DelphMET_Phi)));
        double MET  = tree->Get(MET,"MET");
        double DMET  = tree->Get(DMET,"DelphMET");
        double MT2W  = tree->Get(MT2W,"MT2W");
        double DMT2W  = tree->Get(DMT2W,"MT2W0");
        double MT  = tree->Get(MT,"MT");
        double DMT  = tree->Get(DMT,"MT0");
        
        hbasym->Fill(basym);
        hDphi->Fill(dPhi,EvWeight);
        hDelphDphi->Fill(DelphdPhi,EvWeight);
        hMET->Fill(MET,EvWeight);
        hDelphMET->Fill(DMET,EvWeight);
        hMT2W->Fill(MT2W,EvWeight);
        hDelphMT2W->Fill(DMT2W,EvWeight);
        hMT->Fill(MT,EvWeight);
        if(MT!=0)hMT0->Fill(MT,EvWeight);
        hDelphMT->Fill(DMT,EvWeight);
        if(DMT!=0)hDelphMT0->Fill(DMT,EvWeight);
        
        double MET_Phi  = tree->Get(MET_Phi,"MET_Phi");
        double METx,METy;
        if(D){
            METx = cos(DelphMET_Phi)*DMET;
            METy = sin(DelphMET_Phi)*DMET;
        }
        else{
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
        
        double Meff = MET + HT;
        if(D)  Meff = DMET + HT;
        if(D) { if(DMET>0)hMETMeff->Fill(DMET/Meff,EvWeight); } // don't remove the curly bracket!
        else  if(MET>0)hMETMeff->Fill(MET/Meff,EvWeight);
        
        double AW = Meff/sqrt(HT);
        double AW1 = MT2W/sqrt(HT);
        double AW2 = MT2W/sqrt(HT+MET);
        double HTMET = sqrt(HT+MET);
        double sHT = sqrt(HT);
        double sMET = sqrt(MET);
        double AW3 = MT2W/sqrt(MT+MET);
        double AW4 = MT2W/sqrt(MT+HT+MET);
        
        hAW->Fill(AW, EvWeight);
        hAW1->Fill(AW1, EvWeight);
        hAW2->Fill(AW2, EvWeight);
        hAW3->Fill(AW3, EvWeight);
        hAW4->Fill(AW4, EvWeight);
        hDenom->Fill(HTMET, EvWeight);
        hsHT->Fill(sHT, EvWeight);
        hsMET->Fill(sMET, EvWeight);
        hCentrality->Fill(C,EvWeight);
        
        if(D){
            if(DMET<250) continue;
            hMET_cut->Fill(DMET,EvWeight);
        }
        else{
            if(MET<250) continue;
            hMET_cut->Fill(MET,EvWeight);
        }

        CFCounter[4]+= EvWeight;
        iCFCounter[4]++;
        
        if(D){
            if(DelphdPhi < 0.5) continue;
            hDphi_cut->Fill(DelphdPhi,EvWeight);
        }
        else{
            if(dPhi < 0.5) continue;
            hDphi_cut->Fill(dPhi,EvWeight);
        }
        
        CFCounter[5]+= EvWeight;
        iCFCounter[5]++;
        
        hLepPt_2->Fill(lepPt,EvWeight);
        hLepEta_2->Fill(lepEta,EvWeight);
        
        // 5. MT cut
        if(D){
            if(DMT<120) continue;
            hMT_cut->Fill(DMT,EvWeight);
        }
        else{
            if(MT<120) continue;
            hMT_cut->Fill(MT,EvWeight);
        }
        
        CFCounter[6]+= EvWeight;
        iCFCounter[6]++;
        
        if(D){ if(DMT2W<170) continue; hMT2Wpre->Fill(DMT2W,EvWeight);}
        else if(MT2W<170) continue; hMT2Wpre->Fill(MT2W,EvWeight);
        
        CFCounter[7]+= EvWeight;
        iCFCounter[7]++;
        
        if(HT < 500) continue;
        hHT_cut->Fill(HT,EvWeight);
        
        CFCounter[8]+= EvWeight;
        iCFCounter[8]++;
        
        hLepPt_3->Fill(lepPt,EvWeight);
        hLepEta_3->Fill(lepEta,EvWeight);
        
        HT_CUT.push_back(HT);
        MET_CUT.push_back(MET);
        MT2W_CUT.push_back(MT2W);
        MT_CUT.push_back(MT);
        VEC_AWH2.push_back(AW2);
        VEC_AWH1.push_back(AW1);
        
    }
    
    //TreeS->Fill();
    //TreeB0->Fill();
    //TreeB1->Fill();
    //TreeB2->Fill();
    TreeB3->Fill();
    
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
    
    //Write
    hHT_0->Write();
    hNel->Write();
    hNmu->Write();
    hNjet->Write();
    hLepPt_0->Write();
    hLepEta_0->Write();
    hNjet_check->Write();
    hmaxbJetPt->Write();
    hmaxpbJetPt->Write();
    hNbjet->Write();
    hHT_1->Write();
    hLepPt_1->Write();
    hLepEta_1->Write();
    hJetPt_1->Write();
    hbasym->Write();
    hDphi->Write();
    hDelphDphi->Write();
    hMET->Write();
    hMT2W->Write();
    hDelphMT2W->Write();
    hMT->Write();
    hMT0->Write();
    hDelphMT->Write();
    hDelphMT0->Write();
    hCentrality->Write();
    hMETMeff->Write();
    hAW->Write();
    hAW1->Write();
    hAW2->Write();
    hAW3->Write();
    hAW4->Write();
    hDenom->Write();
    hsHT->Write();
    hsMET->Write();
    hMET_cut->Write();
    hDphi_cut->Write();
    hLepPt_2->Write();
    hLepEta_2->Write();
    hMT_cut->Write();
    hMT2Wpre->Write();
    hHT_cut->Write();
    hLepPt_3->Write();
    hLepEta_3->Write();
    
    TreeS->Write();
    TreeB0->Write();
    TreeB1->Write();
    TreeB2->Write();
    TreeB3->Write();
}
    
    

