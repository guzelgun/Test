#! /usr/bin/env python

from ROOT import *
import sys
from DataFormats.FWLite import Events, Handle
from glob import glob


# Make VarParsing object 
# https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideAboutPythonConfigFile#VarParsing_Example
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')
options.parseArguments()

def main():


    inputFiles = options.inputFiles
    if inputFiles == []:
        print 'running on the EDM file AODSIM'
        inputFiles = ['/nfs/dust/cms/user/guezelgc/disappearingTracks/CMSSW_8_0_21/src/AODSIM/AODSIM10CM/pMSSM12_MCMC1_27_200970_step3_AODSIM_10_9nFiles1.root']


    events = Events(inputFiles)

# create handle and labels of loop

    handle_genparticles  = Handle ("vector<reco::GenParticle>") #strings to use here are learned from AODSIM by edmdumpeventcontent
    label_genparticles = ('genParticlePlusGeant')

    handle_jets = Handle ("vector<reco::PFJet>")
    label_jets = ('ak4PFJetsCHS')

    handle_muons = Handle ("vector<reco::Muon>")
    label_muons = ('muons')

    handle_pfmet = Handle ("vector<reco::PFMET>")
    label_pfmet = ('pfMet')

    handle_tracks  = Handle ("vector<reco::Track>")
    label_tracks = ('generalTracks')

# Create histograms, etc.
    gROOT.SetBatch()        # don't pop up canvases
    gROOT.SetStyle('Plain') # white background 

#define histograms
    hPt                 = TH1F("hPt", "gen #chi^{#pm} Pt", 25, 0, 1000)  # fill pt of gen Chargino
    hPtg                = TH1F("hPtg", "gen Gluino Pt", 25, 0, 1000)  # fill pt of gen Gluino
    hPtgSM              = TH1F("hPtgSM", "gen Gluon Pt", 25, 0, 1000)  # fill pt of gen Gluino
    hGluinoMass		= TH1F("hGluinoMass", "gen Gluino Mass", 25, 0, 2000)
    hCharginoMass	= TH1F("hCharginoMass", "gen #chi^{#pm} Mass ", 25, 0, 1000)
    hHT                 = TH1F("hHT", "HT", 30, 1000, 2000)                 # fill HT of events
    hJets               = TH1F("hJets", "Number of Jets", 20, 0, 20)   # fill number of jets
    hPtMU		= TH1F("hPtMU", "Muon Pt", 25, 0, 150)		 # fill pt of Muons
    hMETJetDR		= TH1F("hMETJetDR", "(#Delta #Phi(MET, Max PT Jet))", 70, -3, 3)		 # fill pt of Muons
    hDrChiTrack         = TH1F("hDrChiTrack", "min #DeltaR(#chi^{#pm},track)", 200, -.001, .01) # fill min Delta R of charginos
    hjetMax		= TH1F("hjetMax", "jet with the highest pt", 200, 1, 200) # fill min Delta R of charginos
    hDecayLength	= TH1F("hDecayLength", "All gen #chi^{#pm} Decay Length", 30, -5, 5)
    hgenLife		= TH1F("hgenLife", "Generation point of Chargino", 30, -1, 3)
    hPfMET		= TH1F("hPfMET", "MET Pt", 25, 0, 1000)  # fill pt of MET

    hdispTrackPt	= TH1F("hdispTrackPt", "pt of the Disappearing Tracks", 200, 1, 500)
    hdispTrackEta	= TH1F("hdispTrackEta", "#eta of the Disappearing Tracks", 70, -3, 3)
    hdispTrackDPhi	= TH1F("hdispTrackDPhi", "(#Delta#Phi(MET, Disappearing Tracks))", 70, -3, 3)

    hbckTrackPt		= TH1F("hbckTrackPt", "pt of the Background Tracks", 200, 1, 500)
    hbckTrackEta	= TH1F("hbckTrackEta", "#eta of the Background Tracks", 70, -3, 3)
    hbckTrackDPhi	= TH1F("hbckTrackDPhi", "(#Delta#Phi(MET, Background Tracks))", 70, -3, 3)

    hValidHits          = TH1F("hValidHits","number of valid hits^{#DeltaR(#chi^{#pm},track) < 0.001}", 15, 0 ,30)
    hInnerMiss          = TH1F("hInnerMiss","missing inner hits^{#DeltaR(#chi^{#pm},track) < 0.001}", 15, 0 ,15)
    hMiddleMiss         = TH1F("hMiddleMiss","missing middle hits^{#DeltaR(#chi^{#pm},track)< 0.001}", 15, 0, 15)
    hOuterMiss          = TH1F("hOuterMiss","missing outer hits^{#DeltaR(#chi^{#pm},track) < 0.001}",15, 0, 15)

    hValidHits2          = TH1F("hValidHits2","number of valid hits^{#DeltaR(#chi^{#pm},track) > 0.002}", 15, 0 ,30)
    hInnerMiss2          = TH1F("hInnerMiss2","missing inner hits^{#DeltaR(#chi^{#pm},track) > 0.002}", 15, 0 ,15)
    hMiddleMiss2         = TH1F("hMiddleMiss2","missing middle hits^{#DeltaR(#chi^{#pm},track) > 0.002}", 15, 0, 15)
    hOuterMiss2          = TH1F("hOuterMiss2","missing outer hits^{#DeltaR(#chi^{#pm},track) > 0.002}",15, 0, 15)
    hDecayLengthTag      = TH1F("hDecayLengthTag","gen #chi^{#pm} with RECO Tracks Decay Length", 30, -5, 5)

   
	

#loop over events
    for event in events:

        event.getByLabel (label_genparticles, handle_genparticles)
        event.getByLabel (label_jets, handle_jets)
	event.getByLabel (label_muons, handle_muons)	
	event.getByLabel (label_pfmet, handle_pfmet)	
        event.getByLabel (label_tracks, handle_tracks)

    # get the product 

        genparticles = handle_genparticles.product()
        jets = handle_jets.product()
	muons = handle_muons.product()
        pfmet = handle_pfmet.product()
        tracks = handle_tracks.product()


#loop over the gen particles, can also loop over muons, electrons jets etc
        for gp in genparticles:

            if gp.pt() < 50 :continue  #uncomment to add pt threshold. Can also add phi, eta threshold.
            if abs(gp.pdgId())==1000024:
                hPt.Fill(gp.pt())          # can look for any other particles via pdgID like gluinos or your favourite particle 

	for gp in genparticles:
            if abs(gp.pdgId())==1000021:
                hPtg.Fill(gp.pt())          # can look for any other particles via pdgID like gluinos or your favourite particle
		hGluinoMass.Fill(gp.mass())
	for gp in genparticles:
            if abs(gp.pdgId())==21:
                hPtgSM.Fill(gp.pt())

        for mu in muons:
                hPtMU.Fill(mu.pt())

	#ptMax=0
	#for jet in jets:
	#    if jet.pt() > ptMax:
	#        ptMax = jet.pt()
	#jetMax.Fill(ptMax)

	ptmax = 0
        for ijet, jet in enumerate(jets):
            if jet.pt() > 0:
                ptMax = jet.pt()
                jet_id = ijet
		

	#AKSHANSH' SCRIPT
        charginos = []
        nChi = 0
        nDT  = 0
        flag = 0 # flag 1 later if DT found
        flagChi = 0
        for gp in genparticles:
            if gp.pt() < 30 :continue # or abs(gp.eta()) > 1.0):continue
            #if abs(gp.pdgId())==1000021 and gp.status()==1:
                #hGluinoMass.Fill(gp.mass())
            if abs(gp.pdgId())==1000024 and gp.status()==1:
                hCharginoMass.Fill(gp.mass())
		#            print 'found chargino!, pt=', gp.pt() , 'mass =' , gp.m()
                #            exit(0)
                flagChi = 1
                trkM  = 0 # number of tracks in dr < .1 / .2/ .3
                trkM2 = 0
                trkM3 = 0
                try:
  	        	log10decaylength = (TMath.Sqrt(pow(gp.daughter(0).vx() - gp.vx(),2) + pow(gp.daughter(0).vy()-gp.vy(),2))) 
                	loglog = TMath.Log10(log10decaylength)
	       		#mLife.Fill(loglog)
		 	log10Genpoint = (TMath.Sqrt(pow(gp.vx(),2) + pow(gp.vy(),2) + pow(gp.vz(),2) )) 
			#log10Genpoint = (TMath.Sqrt(pow(gp.vx(),2) + pow(gp.vy(),2) ))                	
			loglog1 = TMath.Log10(log10Genpoint)
			#genLife.Fill(loglog1)
	    	except:	
			loglog = 10 
          
		hDecayLength.Fill(loglog)
                hPt.Fill(gp.pt())
                #hchieta.Fill(gp.eta())
                #hlog10ctau.Fill(log10decaylength)
                #hpfMET.Fill(MET[0].pt())
                nChi = nChi +1
                chiTlv = TLorentzVector()
                chiTlv.SetPxPyPzE(gp.px(),gp.py(),gp.pz(),gp.energy())
                charginos.append(chiTlv)
                DPtrelmax = 2
                drsmall = .02
                Dptmax  = 3
                track_id = -1
                track_idr = -1
                flagtrk = 0

		#1ST SET HT CUT

		HT = getHT(jets)		

		#if HT < 300: continue
  
		#if (pfmet[0] < 100) or HT <300: continue
		
                  
                for itrack, track in enumerate(tracks):
                    if track.pt() < 15: continue
                    if flagtrk == 0:hDecayLength.Fill(loglog) # fill the gen chargino hist only if atleast one track with Pt > 15 is found
                    flagtrk = 1
                    trkTlv = TLorentzVector()
                    trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt())
                    dr = trkTlv.DeltaR(chiTlv)

                    if dr > 0.002:
		    	hbckTrackPt.Fill(tracks[track_id].pt())
		    	hbckTrackEta.Fill(tracks[track_id].eta())
		    	hbckTrackDPhi.Fill(pfmet[0].phi()-tracks[track_id].phi())
                    	hValidHits2.Fill(tracks[track_id].numberOfValidHits())
                    	hInnerMiss2.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS))
                    	hMiddleMiss2.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS))
                    	hOuterMiss2.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))

                    Dpt = (gp.pt()-track.pt())/(.5*(gp.pt()+track.pt()))
                    DPtrel = DPtRel(track,gp)   # get relative err on Pt
                    if dr < .1 : trkM  = trkM  + 1
                    if dr < .2 : trkM2 = trkM2 + 1
                    if dr < .3 : trkM3 = trkM3 + 1
                    if dr<drsmall:
                        drsmall = dr
                        track_id = itrack
                        DptDRmatch = (gp.pt()-track.pt())/(.5*(gp.pt()+track.pt()))
                    if drsmall < 0.02:
                        if abs(Dpt)< abs(Dptmax):
                            Dptmax = Dpt
                            track_idd = itrack
                        if DPtrel < DPtrelmax:
                            DPtrelmax = DPtrel
                            track_idr = itrack
                                    
                hDrChiTrack.Fill(drsmall)
		if drsmall < 0.006:
		    hdispTrackPt.Fill(tracks[track_id].pt())
		    hdispTrackEta.Fill(tracks[track_id].eta())
		    hdispTrackDPhi.Fill(pfmet[0].phi()-tracks[track_id].phi())
		
                    hValidHits.Fill(tracks[track_id].numberOfValidHits())
                    hInnerMiss.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS))
                    hMiddleMiss.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS))
                    hOuterMiss.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))
		    hDecayLengthTag.Fill(loglog)		     
	#AKSHANSH' SCRIPT        

        #HT = getHT(jets)    # get HT from function  getHT , pass the array of jets. 
        hHT.Fill(HT)
        hPfMET.Fill(pfmet[0].pt())
	hMETJetDR.Fill(pfmet[0].phi()-jets[jet_id].phi())
	nJets = getJET(jets)
        hJets.Fill(nJets)

    #fnew = TFile('sample_DENEME.root','recreate') # save hists in a root file
    identifier = inputFiles[0][inputFiles[0].rfind('/')+1:].replace('.root','').replace('_step2','').replace('_AODSIM','').replace('_*','').replace('*','')
    identifier+='nFiles'+str(len(inputFiles))
    fnew = TFile('hists_'+identifier+'.root','recreate')

    hPt.Write()
    hHT.Write()
    hJets.Write()
    hPtMU.Write()
    hDrChiTrack.Write()
    hPtg.Write()
    hPtgSM.Write()
    hPfMET.Write()
    hMETJetDR.Write()
    hdispTrackPt.Write()
    hdispTrackEta.Write()
    hdispTrackDPhi.Write()
    hValidHits.Write()
    hInnerMiss.Write()
    hMiddleMiss.Write()
    hOuterMiss.Write()
    hbckTrackPt.Write()
    hbckTrackEta.Write()
    hbckTrackDPhi.Write()
    hValidHits2.Write()
    hInnerMiss2.Write()
    hMiddleMiss2.Write()
    hOuterMiss2.Write()
    hDecayLength.Write()
    hgenLife.Write()
    hCharginoMass.Write()
    hGluinoMass.Write()
    hDecayLengthTag.Write()

#Save the histograms as pdf or png
    
    c1 = TCanvas()
    hPt.Draw()
    #c1.Print("hPt.pdf")

    hHT.Draw()
    #c1.Print("hHT.pdf")

    hPtg.Draw()
    #c1.Print("hPtg.pdf")

    hPtgSM.Draw()
    #c1.Print("hPtgSM.pdf")

    hJets.Draw()
    #c1.Print("nJets.pdf")

    hPtMU.Draw()
    #c1.Print("hPtMU.pdf")

    hPfMET.Draw()
    #c1.Print("hPfMET1.pdf")

    hDrChiTrack.Draw()
    #c1.Print("hDrChiTrack.pdf")

    #jetMax.Draw()
    #c1.Print("jetMax.pdf")

def DPtRel(track, genp):
    dptrel = abs((track.pt()-genp.pt())/genp.pt())
    return dptrel

def getHT(jetlist = [], *args):
    ht =0
    for jet in jetlist:
        if jet.pt()< 25:continue
        ht = ht + jet.pt()
    return ht

def getJET(jetlist = [], *args):
    j =0
    for jet in jetlist:
        if jet.pt() < 20:continue    
        looseID= 1
        if abs(jet.eta()) < 2.4: looseID = jet.neutralEmEnergyFraction() < .99 and jet.neutralHadronEnergyFraction() < .99 and jet.numberOfDaughters() > 1 and jet.chargedEmEnergyFraction() < .99 and jet.chargedHadronEnergyFraction() > 0 and  jet.chargedMultiplicity() > 0
        if abs(jet.eta()) >= 2.4: looseID = jet.neutralEmEnergyFraction() < .99 and jet.neutralHadronEnergyFraction() < .99 and jet.numberOfDaughters() > 1
        if looseID == 0:continue
        j = j + 1    
    return j

main()



