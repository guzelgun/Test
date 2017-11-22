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
    # Events takes either
    # - single file name
    # - list of file names
    # - VarParsing options
    # use Varparsing object
#events = Events (options)
    inputFiles = options.inputFiles
    if inputFiles == []:
        print 'running on the default 100 events sample cTAu= 55 cm'
        inputFiles = ['/nfs/dust/cms/user/guezelgc/disappearingTracks/pMSSM12_MCMC1_27_200970_step2_AODSIM.root']
#        inputFiles = ['/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_10_374794_step2_AODSIM.root']
    events = Events(inputFiles)
#events = Events('/nfs/dust/cms/user/beinsam/LongLiveTheChi/aodsim/smallchunks/pMSSM12_MCMC1_10_374794_step2_AODSIM_2of100.root')
#events = Events('root://cmsxrootd.fnal.gov//store/user/sbein/LongLiveTheChi/pMSSM12_MCMC1_27_200970_step2_AODSIM_15of100.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_10_374794_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_12_865833_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_13_547677_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_20_690321_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_22_237840_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_24_345416_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_27_200970_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_27_969542_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_28_737434_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_37_569964_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_4_252033_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_44_855871_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_47_872207_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_5_448429_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/merged/pMSSM12_MCMC1_8_373637_step2_AODSIM.root')
#events = Events('/nfs/dust/cms/user/singha/MET_scan_8/CMSSW_8_0_20/src/MetScanning/LLPtupleTEST1/AODI/pMSSM12_MCMC1_15sp_step2_AODSIM.root')

#events = glob('LLPtupleTEST1/AODI/merged/pMSSM*.root')

    handle_muons  = Handle ("std::vector<reco::Muon>")
    label_muons = ('muons')

    handle_tracks  = Handle ("vector<reco::Track>")
    label_tracks = ('generalTracks')

    handle_pfcands  = Handle ("std::vector<reco::PFCandidate>")
    label_pfcands = ('particleFlow')

    handle_genparticles  = Handle ("vector<reco::GenParticle>")
    label_genparticles = ('genParticlePlusGeant')

    handle_jets = Handle ("vector<reco::PFJet>")
    label_jets = ('ak4PFJetsCHS')

    handle_MET = Handle ("vector<reco::PFMET>")
    label_MET = ('pfMet')


#print dir(generalTracks)                                                                                                                                                      
#hitpattern = track.hitPattern()                                                                                                                                       
#print dir(hitpattern)                                                                                                                                                 
#exit(0)


# Create histograms, etc.
    gROOT.SetBatch()        # don't pop up canvases
    gROOT.SetStyle('Plain') # white background
    gStyle.SetOptStat(111111)

    hchieta             = TH1F("hchieta", "#chi #eta distribution", 100, -3.5,3.5 )
    hDrChiTrack         = TH1F("hDrChiTrack", "min #DeltaR(#chi^{#pm},track)", 50, -0.00001, .5)
    hDrRandomTrackTrack = TH1F("hDrRandomTrackTrack", "hDrRandomTrackTrack", 100, -0.00001, 1)
    
    hdptrelChiTrack     = TH1F("hdptrelChiTrack", "DPtRel(#chi^{#pm},track)", 120, 0, 1.2) 


    hDptChiTrack        = TH1F("hDptChiTrack", "min DPt(#chi^{#pm},track)^{#DeltaR(#chi^{#pm},track)<0.02}", 200, -2, 2)

    hDptChiTrackDRmatch = TH1F("hDptChiTrackDRmatch", "DPt(#chi^{#pm},track)^{min #DeltaR(#chi^{#pm},track)}", 200, -2, 2)
    hpfgen              = TH1F("hpfgen", "Wrongly matched SM particles", 10, 0, 10)

    hlog10ctau          = TH1F("hlog10ctau", "#chi^{#pm} #DeltaS (gen level)", 100,.1,1000.0)
    hlog10ctauTag       = TH1F("hlog10ctauTag", "#chi^{#pm} #DeltaS(tagged tracks)", 100,.1,1000.0)

    hpfMET              = TH1F("hpfMET", "PFMet", 30, 0, 1500)
    hpfMETTag           = TH1F("hpfMETTag", "PFMet(events with disappearing tracks)", 30, 0, 1500)

    hpfMETsample        = TH1F("hpfMETsample", "PFMet(events simulated)", 30, 0, 1500)
    hpfMETsampleTag     = TH1F("hpfMETsampleTag", "PFMet(events with charginos)", 30, 0, 1500)


## Track hit pattern variables
    hnDT                = TH1F("hnDT", "number of disappearing tracks", 4, 0, 4)
    hnChi               = TH1F("hnChi", "number of charginos", 4, 0, 4)

    hValidHits          = TH1F("hValidHits","number of valid hits^{#DeltaR(#chi^{#pm},track)<0.02}", 15, 0 ,30)
    hInnerMiss          = TH1F("hInnerMiss","missing inner hits^{#DeltaR(#chi^{#pm},track)<0.02}", 15, 0 ,15)
    hMiddleMiss         = TH1F("hMiddleMiss","missing middle hits^{#DeltaR(#chi^{#pm},track)<0.02}", 15, 0, 15)
    hOuterMiss          = TH1F("hOuterMiss","missing outer hits^{#DeltaR(#chi^{#pm},track)<0.02}",15, 0, 15)
    
    hValidHitsDT        = TH1F("hValidHitsDT","number of valid hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0 ,30)
    hInnerMissDT        = TH1F("hInnerMissDT","missing inner hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0 ,15)
    hMiddleMissDT       = TH1F("hMiddleMissDT","missing middle hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0, 15)
    hOuterMissDT        = TH1F("hOuterMissDT","missing outer hits^{#DeltaR(#chi^{#pm},track)<0.01}",15, 0, 15)
    
##More track hit pattern Variables
    htrkHits            = TH1F("htrkHits","number of tracker hits^{#DeltaR(#chi^{#pm},track)<0.01}", 30, 0 ,30)
    hpixHits            = TH1F("hpixHits","number of pixel hits^{#DeltaR(#chi^{#pm},track)<0.01}", 8, 0 ,8)
    hstripHits          = TH1F("hstripHits","number of strip hits^{#DeltaR(#chi^{#pm},track)<0.01}", 20, 0 ,20)
    hpixBHits           = TH1F("hpixBHits","number of pixel barrel hits^{#DeltaR(#chi^{#pm},track)<0.01}", 5, 0 ,5)
    hpixECHits          = TH1F("hpixECHits","number of pixel end cap hits^{#DeltaR(#chi^{#pm},track)<0.01}", 4, 0 ,4)
    hstripTIBhits       = TH1F("hstripTIBhits","number of inner barrel hits^{#DeltaR(#chi^{#pm},track)<0.01}", 10, 0 ,10)
    hstripTIDhits       = TH1F("hstripTIDhits","number of inner disk hits^{#DeltaR(#chi^{#pm},track)<0.01}", 10, 0 ,10)
    hstripTOBhits       = TH1F("hstripTOBhits","number of outer barrel hits^{#DeltaR(#chi^{#pm},track)<0.01}", 10, 0 ,10)
    hstripTEChits       = TH1F("hstripTEChits","number of tracker end cap hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0 ,15)
    hOutmissBar         = TH1F("hOutmissBar","barrel(#eta#leq.9) missing outer hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0 ,15)
    hOutmissMix         = TH1F("hOutmissMix",".9 <#eta #leq2.1 missing outer hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0 ,15)
    hOutmissEC          = TH1F("hOutmissEC","end cap(#eta #geq 2.1) missing outer hits^{#DeltaR(#chi^{#pm},track)<0.01}", 15, 0 ,15)
### event selection vars
    hdphi               = TH1F("hdphi", "#Delta #phi (#chi,jet^{1})", 30, 0, 3.15)
    hHT                 = TH1F("hHT", "HT", 30, 0, 1500)
    hHTTag              = TH1F("hHTTag", "HT (events with disapperaing tracks)", 30, 0, 1500)
    hnbjet              = TH1F("hnbjet", "b jet multiplicity", 10, 0, 10)
    hnbjetTag           = TH1F("hnbjetTag", "b jet multiplicity (events with disapperaing tracks)", 10, 0, 10)
    hnjet               = TH1F("hnjet", "jet multiplicity", 10, 0, 12)
    hnjetTag            = TH1F("hnjetTag", "jet multiplicity (events with disapperaing tracks)", 10, 0, 12)
    hPt                 = TH1F("hPt", "gen #chi^{#pm} Pt", 25, 0, 1000)
    hPtTag              = TH1F("hPtTag", "tagged tracks #chi^{#pm} Pt", 25, 0, 1000)
###Isolation
    htrkdR              = TH1F("htrkdR", "number of tracks in #DeltaR < .01 ", 10, 0, 10)
    htrkdR2             = TH1F("htrkdR2", "number of tracks in #DeltaR < .02 ", 10, 0, 10)
    hrelisop3           = TH1F("hrelisop3", "RelIso #DeltaR < 0.3 ", 30, -0.0001, 0.5)

#    print 'track pt', 'gen particle pt', 'gen pdgId', 'dR(trk,gen)','dr(trk,chi)'
# loop over events
    for event in events:
        # use getByLabel, just like in cmsRun
        event.getByLabel (label_muons, handle_muons)
        event.getByLabel (label_tracks, handle_tracks)
        event.getByLabel (label_genparticles, handle_genparticles)
        event.getByLabel (label_pfcands, handle_pfcands)
        
        event.getByLabel (label_MET, handle_MET)
        event.getByLabel (label_jets, handle_jets)
    # get the product
        muons = handle_muons.product()
        tracks = handle_tracks.product()
        genparticles = handle_genparticles.product()
        pfcands = handle_pfcands.product()
        
        MET = handle_MET.product()
        jets = handle_jets.product()
        
        #    print dir(MET)
        #    print '========'
        #    print dir(jets)
        #    print '========'
        
        charginos = []
        #    print '='*10
        #    print dir(tracks)
        #    hitpattern = tracks.hitPattern() 
        #    print dir(hitpattern)
        #    exit(0)
        
#        print 'check 2'
        
#    print dir(MET[0])
        
    #    print MET[0].pt()
        nChi = 0
        nDT  = 0
        hpfMETsample.Fill(MET[0].pt())
        flag = 0 # flag 1 later if DT found
        flagChi = 0
        for gp in genparticles:
            if gp.pt() < 30 :continue # or abs(gp.eta()) > 1.0):continue
            if abs(gp.pdgId())==1000024 and gp.status()==1:
                #            print 'found chargino!, pt=', gp.pt() , 'mass =' , gp.m()
                #            exit(0)
                flagChi = 1
                trkM  = 0 # number of tracks in dr < .1 / .2/ .3
                trkM2 = 0
                trkM3 = 0
                try:
                    log10decaylength = (TMath.Sqrt(pow(gp.daughter(0).vx() - gp.vx(),2) + pow(gp.daughter(0).vy()-gp.vy(),2) )) #TMath.Log10

                #print 'Daughter 1 Pt:', gp.daughter(0).pt(), 'pdgID:', gp.daughter(0).pdgId()
                #print 'Daughter 2 Pt:', gp.daughter(1).pt(), 'pdgID:', gp.daughter(1).pdgId()
                #print 'Chargino Pt  :', gp.pt(), 'pdgID:', gp.pdgId()
                #print '='*10
                except:
                #print 'no daughters'
                    log10decaylength = -1.49
          
                hPt.Fill(gp.pt())
                hchieta.Fill(gp.eta())
               # hlog10ctau.Fill(log10decaylength)
            #     hpfMET.Fill(MET[0].pt())
                nChi = nChi +1
                chiTlv = TLorentzVector()
                chiTlv.SetPxPyPzE(gp.px(),gp.py(),gp.pz(),gp.energy())
                charginos.append(chiTlv)
                DPtrelmax = 2
                drsmall = .02
                Dptmax     = 3
                track_id = -1
                track_idr = -1
                flagtrk = 0
                for itrack, track in enumerate(tracks):
                    if track.pt() < 15: continue
                    if flagtrk == 0:hlog10ctau.Fill(log10decaylength) # fill the gen chargino hist only if atleast one track with Pt > 15 is found
                    flagtrk = 1
                    trkTlv = TLorentzVector()
                    trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt())
                    dr = trkTlv.DeltaR(chiTlv)
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
                                    
                hDrChiTrack.Fill(drsmall)    # min dr between gen Chi and tracks, one per chargino
                hDptChiTrack.Fill(Dptmax)    # min Dpt between gen Chi and tracks(only the tracks with dr < .02), one per event
                hdptrelChiTrack.Fill(DPtrel) # 
                                #            print 'dPt', DptDRmatch
                                #            print 'dR', drsmall
                                #            print 'BM', tracks[track_idr].pt()
                                #            print 'BM', tracks[track_id].pt()
                                #            print 'chi pt', gp.pt()
                                
            #           print '='*10
                htrkdR.Fill(trkM)
                htrkdR2.Fill(trkM2)
                if drsmall < 0.02:
                        
                    relErr = tracks[track_id].ptError()/tracks[track_id].pt()

#=================
#trackerLayersWithoutMeasurement
#numberOfHits
                    hDptChiTrackDRmatch.Fill(DptDRmatch) # Dpt of best dR match

                    hValidHits.Fill(tracks[track_id].numberOfValidHits())
                    hInnerMiss.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS))
                    hMiddleMiss.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS))
                    hOuterMiss.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))
                    
                if drsmall < 0.01: 
                    riso = reliso(tracks[track_id],tracks)# and abs(Dptmax) < 1.0:

                if drsmall < 0.01:# and riso < 0.1:
#                    riso = reliso(tracks[track_id],tracks)
#                    hdSvsIso.Fill(riso,log10decaylength)
                    hrelisop3.Fill(riso)
#                    print reliso(tracks[track_id],tracks)
                    hlog10ctauTag.Fill(log10decaylength)
                    hValidHitsDT.Fill(tracks[track_id].numberOfValidHits())
                    hInnerMissDT.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS))
                    hMiddleMissDT.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS))
                    hOuterMissDT.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))

######More Track Vars
                    htrkHits.Fill(tracks[track_id].hitPattern().trackerLayersWithMeasurement())
                    hpixHits.Fill(tracks[track_id].hitPattern().pixelLayersWithMeasurement())
                    hstripHits.Fill(tracks[track_id].hitPattern().stripLayersWithMeasurement())
                    hpixBHits.Fill(tracks[track_id].hitPattern().pixelBarrelLayersWithMeasurement())
                    hpixECHits.Fill(tracks[track_id].hitPattern().pixelEndcapLayersWithMeasurement())
                    hstripTIBhits.Fill(tracks[track_id].hitPattern().stripTIBLayersWithMeasurement())
                    hstripTIDhits.Fill(tracks[track_id].hitPattern().stripTIDLayersWithMeasurement())
                    hstripTOBhits.Fill(tracks[track_id].hitPattern().stripTOBLayersWithMeasurement())
                    hstripTEChits.Fill(tracks[track_id].hitPattern().stripTECLayersWithMeasurement())
                
                    if abs(tracks[track_id].eta()) <= .9:
                        hOutmissBar.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))
                    elif abs(tracks[track_id].eta()) > .9 and abs(tracks[track_id].eta()) < 2.1:
                        hOutmissMix.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))
                    elif abs(tracks[track_id].eta()) >= 2.1:
                        hOutmissEC.Fill(tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS))
             
                    flag = 1
                    nDT = nDT + 1
                    dphi = oldFriend(tracks[track_id], jets[0])
                    hdphi.Fill(dphi)
                    hPtTag.Fill(gp.pt())
                    for gen in genparticles:
                        if (gen.pt() < 30 or abs(gen.eta()) > 1.0):continue
                        genTlv = TLorentzVector()
                        genTlv.SetPxPyPzE(gen.px(),gen.py(),gen.pz(),gen.energy())
                        tkTlv  = TLorentzVector()
                        tkTlv.SetPxPyPzE(tracks[track_id].px(),tracks[track_id].py(),tracks[track_id].pz(),tracks[track_id].pt())
                        DR = tkTlv.DeltaR(genTlv)
                        if DR < 0.01 :
                            if abs(gen.pdgId()) == 1000024:hpfgen.Fill(9)
                            if abs(gen.pdgId()) == 211:hpfgen.Fill(3)
                            elif abs(gen.pdgId()) < 3000:hpfgen.Fill(5)
                            #print tracks[track_id].pt(),'*', gen.pt(),'*', gen.pdgId(),'*', DR, '*',drsmall 

        if flagChi ==1:#flagChi
            hnChi.Fill(nChi)
            hnDT.Fill(nDT)
            hpfMETsampleTag.Fill(MET[0].pt())
            hpfMET.Fill(MET[0].pt())
            HT = getHT(jets)
            hHT.Fill(HT)
            njet = getnjets(jets)
            hnjet.Fill(njet)
            nbjet = getnbjets(jets)
            hnbjet.Fill(nbjet)
        if flag ==1:
            hpfMETTag.Fill(MET[0].pt())
#            hnChivsnDT.Fill(nDT,nChi)
            HTtag = getHT(jets)
            hHTTag.Fill(HTtag)
            njettag = getnjets(jets)
            hnjetTag.Fill(njettag)
            nbjettag = getnbjets(jets)
            hnbjetTag.Fill(nbjettag)

           # print(test(gp,tracks,jets))
#####testing 

#       print 'relax, everything is fine :)'

    identifier = inputFiles[0][inputFiles[0].rfind('/')+1:].replace('.root','').replace('_step2','').replace('_AODSIM','').replace('_*','').replace('*','')
    identifier+='nFiles'+str(len(inputFiles))
    
    fnew = TFile('histsEDM_'+identifier+'.root','recreate')
#    fnew = TFile('analiseTEST.root','recreate')

    hValidHits.Write()
    hInnerMiss.Write()
    hMiddleMiss.Write()
    hOuterMiss.Write()

    hValidHitsDT.Write()
    hInnerMissDT.Write()
    hMiddleMissDT.Write()
    hOuterMissDT.Write()

    hDrChiTrack.Write()
#hDrRandomTrackTrack.Write()
#c1.Print ("dr_track_randomtrack.pdf")
    hDptChiTrack.Write()

    hlog10ctau.Write()
    hlog10ctauTag.Write()
######Isolation plots

    hrelisop3.Write()
#    hdSvsIso.Write()
#    hetavsIso.Write()
##event selection vars

    hdphi.Write()
    hnjet.Write()
    hnjetTag.Write()
    hnbjet.Write()
    hnbjetTag.Write()

    hpfMET.Write()
    hpfMETTag.Write()
    hHT.Write()
    hHTTag.Write()
#####MC verification and cut optimization
    htrkdR2.Write()
    htrkdR.Write()
    hpfgen.Write()

    hdptrelChiTrack.Write()
    hpfMETsample.Write()
    hpfMETsampleTag.Write()
    hnDT.Write()
    hnChi.Write()
    hPt.Write()
    hPtTag.Write()

    htrkHits.Write()
    hpixHits.Write()
    hstripHits.Write()

    hpixBHits.Write()
    hpixECHits.Write()

    hstripTIBhits.Write()
    hstripTIDhits.Write()
    hstripTOBhits.Write()
    hstripTEChits.Write()

    hOutmissBar.Write()
    hOutmissMix.Write()
    hOutmissEC.Write()
    hchieta.Write()

def reliso(trk,trklist =[], *args):
    trkTlv = TLorentzVector()
    trkTlv.SetPxPyPzE(trk.px(), trk.py(), trk.pz(), trk.pt())
    trackiso = 0
    sumPt = 0
    for track in trklist:
        trackTlv = TLorentzVector()
        trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt())
        if trackTlv.DeltaR(trkTlv) < 0.3 : sumPt = sumPt + track.pt()
    trackiso = (sumPt - (trk.pt()))/trk.pt()
    return trackiso
        
def getHT(jetlist = [], *args):
    ht =0
    for jet in jetlist:
        if jet.pt()< 25:continue
        looseID = 1
        if abs(jet.eta()) < 2.4: looseID = jet.neutralEmEnergyFraction() < .99 and jet.neutralHadronEnergyFraction() < .99 and jet.numberOfDaughters() > 1 and jet.chargedEmEnergyFraction() < .99 and jet.chargedHadronEnergyFraction() > 0 and  jet.chargedMultiplicity() > 0
        if abs(jet.eta()) >= 2.4: looseID = jet.neutralEmEnergyFraction() < .99 and jet.neutralHadronEnergyFraction() < .99 and jet.numberOfDaughters() > 1
        if looseID == 0:continue

        ht = ht + jet.pt()
    return ht
def getnjets(jetlist = [], *args):
    njets =0
    for jet in jetlist:
        if jet.pt()< 25:continue

        looseID= 1
        if abs(jet.eta()) < 2.4: looseID = jet.neutralEmEnergyFraction() < .99 and jet.neutralHadronEnergyFraction() < .99 and jet.numberOfDaughters() > 1 and jet.chargedEmEnergyFraction() < .99 and jet.chargedHadronEnergyFraction() > 0 and  jet.chargedMultiplicity() > 0
        if abs(jet.eta()) >= 2.4: looseID = jet.neutralEmEnergyFraction() < .99 and jet.neutralHadronEnergyFraction() < .99 and jet.numberOfDaughters() > 1
        if looseID == 0:continue

        njets = njets + 1
    return njets
def getnbjets(jetlist = [], *args):
    nbjets =0
    for jet in jetlist:
        if jet.pt()< 25:continue
#        print jet.pdgId()
        if jet.pdgId() == 5:nbjets = nbjets + 1
    return nbjets

def oldFriend(track, jet):
    delPhi = abs (track.phi()- jet.phi())
#    print dir(jet)
#    exit(0)
    return delPhi

def DPtRel(track, genp):
    dptrel = abs((track.pt()-genp.pt())/genp.pt())
    return dptrel

def test(genp,trk =[], jt =[], *args):
#    print 'testing'
#    print trk[0].pt()
#    print jt[0].pt()
#    print genp.pt()
    return 0
main()

   # print dir(track)
   # hitpattern = track.hitPattern()
   # print dir(hitpattern)
   # exit(0)

'''

     #########pixel layer = pix barrel + pixel endcap: test passed                                                                                                              

                print 'tracker layers:',tracks[track_id].hitPattern().trackerLayersWithMeasurement()
                print 'strip layers :', tracks[track_id].hitPattern().stripLayersWithMeasurement()
                l = tracks[track_id].hitPattern().stripTIBLayersWithMeasurement() + tracks[track_id].hitPattern().stripTIDLayersWithMeasurement() + tracks[track_id].hitPattern\
().stripTOBLayersWithMeasurement() + tracks[track_id].hitPattern().stripTECLayersWithMeasurement() + tracks[track_id].hitPattern().pixelLayersWithMeasurement()
                print 'Sum pixel + barrel+ID+EC :', l
                print 'tracker layers without measurement:', tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS), tracks[track_id\
].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS), tracks[track_id].hitPattern().trackerLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'strip without measurement', tracks[track_id].hitPattern().stripLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS), tracks[track_id].hitPattern\
().stripLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS), tracks[track_id].hitPattern().stripLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'pixel without measurement', tracks[track_id].hitPattern().pixelLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS), tracks[track_id].hitPattern\
().pixelLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS), tracks[track_id].hitPattern().pixelLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'TIB without measurement', tracks[track_id].hitPattern().stripTIBLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS) , tracks[track_id].hitPatte\
rn().stripTIBLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS) , tracks[track_id].hitPattern().stripTIBLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'TOB without measurement', tracks[track_id].hitPattern().stripTOBLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS) , tracks[track_id].hitPatte\
rn().stripTOBLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS) , tracks[track_id].hitPattern().stripTOBLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'TID without measurement', tracks[track_id].hitPattern().stripTIDLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS) , tracks[track_id].hitPatte\
rn().stripTIDLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS) , tracks[track_id].hitPattern().stripTIDLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'TEC without measurement', tracks[track_id].hitPattern().stripTECLayersWithoutMeasurement(reco.HitPattern.MISSING_INNER_HITS) , tracks[track_id].hitPatte\
rn().stripTECLayersWithoutMeasurement(reco.HitPattern.TRACK_HITS) , tracks[track_id].hitPattern().stripTECLayersWithoutMeasurement(reco.HitPattern.MISSING_OUTER_HITS)
                print 'track pt, eta, phi', tracks[track_id].pt() , tracks[track_id].eta(), tracks[track_id].phi()


                print '='*10
           #     else:exit(0)                                                                                                                                                   


'''
