#!/usr/bin/python


import sys
from ROOT import gROOT

# scenarios and samples
scenarios = ['NoPU','50PU','140PU']
samples   = ['DiBoson','BosonJets','TopJets','TTbar','483_300','566_375','800_450','1100_450','1125_25']
# by historic reason we have as dir names: bjets  diboson  tdr tdrNew tjets  ttbar for snowmass
# by historic reason we have as dir names: bjets  diboson  tdr tdr tjets  ttbar for snowmass2

#samples w.o weight
noWeight  = ['483_300','566_375','800_450','1100_450','1125_25']

# nTuple input path
base = '/media/cakir/cakir_hdd/snowmass2/'
base2 = '/media/cakir/cakir_hdd/snowmass2/'
# lumi in fb^-1
Lumi=3000# changed for 140PU to 3000


def help():
	print 'First argument analysis:'
	print '   SingleS'
	print ' then a list of inputs - read the code'
	sys.exit(0)

from operator import mul
def scale(fac,list):
	return map(mul,len(list)*[fac],list)

# choose the analysis and a sample
if len(sys.argv)>1:
	if sys.argv[1]=='SingleS':       # single lepton  stop - CMS
		gROOT.LoadMacro('readerSingleS_dmcutflows.C+')
		from ROOT import readerSingleS_dmcutflows as reader
	elif sys.argv[1]=='SingleS_Extract':  # single lepton testing version 
		gROOT.LoadMacro('readerSingleS_Extract.C+')
		from ROOT import readerSingleS_Extract as reader
	elif sys.argv[1]=='SingleS_P':  # single lepton testing version 
		gROOT.LoadMacro('readerSingleS_P.C+')
		from ROOT import readerSingleS_P as reader
	elif sys.argv[1]=='SingleSDelphMET':  # single lepton testing version 
		gROOT.LoadMacro('readerSingleSDelphMET.C+')
		from ROOT import readerSingleSDelphMET as reader
	elif sys.argv[1]=='EWKino':  # EWKino
		gROOT.LoadMacro('readerEWKino.C+') # used to be EWKinoNew
		from ROOT import readerEWKino as reader
	elif sys.argv[1]=='AtlasH':  # Atlas hadronic
		gROOT.LoadMacro('readerAtlasH.C+') # used to be ?????
		from ROOT import readerAtlasH as reader
        elif sys.argv[1]=='Hadronic':  # inclusive hadronic
                gROOT.LoadMacro('readerHadronic.C+') # used to be ?????
                from ROOT import readerHadronic as reader
        elif sys.argv[1]=='Superpol':  # inclusive hadronic
                gROOT.LoadMacro('readerSuperpol.C+') # used to be ?????
                from ROOT import readerSuperpol as reader
        elif sys.argv[1]=='HadronicSUSY':  # inclusive hadronic
                gROOT.LoadMacro('readerHadronicSUSY.C+') # used to be ?????
                from ROOT import readerHadronicSUSY as reader
        else: 
		help()
else: 
	help()


#prepare empty dictionaries
do={}
for scene in scenarios:
	do[scene]={}
	for samp in samples:
		do[scene][samp] = 0
# read input line and set flags what to process
flag=False
for e in sys.argv:
	for scene in scenarios:
		for samp in samples:
			if e==scene+'_'+samp: 
 				do[scene][samp] = 1
 				if scene=='140PU': Lumi=Lumi*10  # lumi*10 for 140PU
				flag=True
if not flag: help()


# prepare empty dictionaries
dirsHT={}
inDir={}
weights={}
for scene in scenarios:
	inDir[scene]={}
#	weights[scene]={}

#--------------------------------------------- sample properties for processing <<<<<<<<<<<<<<<<<<<<<<<<<<<
#
#------------------------------------------------ diboson
# HT dirs 
dirsHT['DiBoson']  = ['0-300/','300-700/','700-1300/','1300-2100/','2100-100000/']
#    xsec
weights['DiBoson'] = [249.97, 35.23, 4.137, 0.417, 0.0477]
#    multiply all by lumi
weights['DiBoson'] = scale(Lumi,weights['DiBoson'])
#
# NoPU
inDir['NoPU']['DiBoson']  = base+'/NoPU/diboson/'
# 50PU
inDir['50PU']['DiBoson']  = base+'/50PU/diboson/'
# 140PU
inDir['140PU']['DiBoson'] = base+'/140PU/diboson/'
#
#------------------------------------------------ boson+jets
# HT dirs 
dirsHT['BosonJets']  = ['0-300/','300-600/','600-1100/','1100-1800/','1800-2700/','2700-3700/','3700-100000/']
#    xsec
weights['BosonJets'] = [34409.92339,2642.85309,294.12311,25.95000,2.42111,0.22690,0.02767]
#    multiply all by lumi
weights['BosonJets'] = scale(Lumi,weights['BosonJets'])
#
# NoPU
inDir['NoPU']['BosonJets']  = base+'/NoPU/bjets/'
# 50PU
inDir['50PU']['BosonJets']  = base+'/50PU/bjets/'
# 140PU
inDir['140PU']['BosonJets'] = base+'/140PU/bjets/'
#
#---------------------------------------------- ttbar
# HT dirs 
dirsHT['TTbar']  = ['0-600/','600-1100/','1100-1700/','1700-2500/','2500-100000/']
#    xsec
weights['TTbar'] = [530.89358,42.55351,4.48209,0.52795,0.05449]
#    multiply all by lumi
weights['TTbar'] = scale(Lumi,weights['TTbar'])
#
# NoPU
inDir['NoPU']['TTbar']  = base+'/NoPU/ttbar/'
# 50PU
inDir['50PU']['TTbar']  = base+'/50PU/ttbar/'
# 140PU
inDir['140PU']['TTbar'] = base+'/140PU/ttbar/'
#
#--------------------------------------------- tjets
# HT dirs 
dirsHT['TopJets']  = ['0-500/','500-1000/','1000-1600/','1600-2400/','2400-100000/']
#    xsec
weights['TopJets'] = [109.73602,5.99325,0.37680,0.03462,0.00312]
#    multiply all by lumi
weights['TopJets'] = scale(Lumi,weights['TopJets'])
#
# NoPU
inDir['NoPU']['TopJets']  = base+'/NoPU/tjets/'
# 50PU
inDir['50PU']['TopJets']  = base+'/50PU/tjets/'
# 140PU
inDir['140PU']['TopJets'] = base+'/140PU/tjets/'
#---------------------------------------------STCs
# HT dirs - none for STC samples
dirsHT['1125_25'] = ['/']
dirsHT['483_300'] = ['/']
dirsHT['566_375'] = ['/']
dirsHT['800_450'] = ['/']
dirsHT['1100_450'] = ['/']
#    xsec - we multiply by lumi

#weights['483_300'] = [2.8]
#weights['566_375'] = [1.2]
#weights['800_450'] = [0.14]
#weights['1100_450'] = [0.014]
#weights['1125_25'] = [0.012048]

weights['483_300'] = [0.8101815]
weights['566_375'] = [0.32587]
weights['800_450'] = [0.0381673]
weights['1100_450'] = [0.00439055]
weights['1125_25'] = [0.00373237]

weights['483_300'] = scale(Lumi,weights['483_300'])
weights['566_375'] = scale(Lumi,weights['566_375'])
weights['800_450'] = scale(Lumi,weights['800_450'])
weights['1100_450'] = scale(Lumi,weights['1100_450'])
weights['1125_25'] = scale(Lumi,weights['1125_25'])

#weights['1125_25'] = [0.012048*Lumi/700000]
#weights['483_300'] = [ 2.8*Lumi/1200000]
#weights['566_375'] = [ 1.2*Lumi/1200000]
#weights['800_450'] = [ 0.14*Lumi/1000000]
#weights['1100_450'] = [ 0.014*Lumi/200000]

#
bit=''
if base2=='/nfs/dust/cms/user/cakir/ECFA_SETUP/SAMPLES_191114/STC8_Delphes_Ntuple/':bit='New' 
# NoPU
inDir['NoPU']['483_300']  = base+'/NoPU/ML/stop_stop/483_300'
inDir['NoPU']['566_375']  = base+'/NoPU/ML/stop_stop/566_375'
inDir['NoPU']['800_450']  = base+'/NoPU/ML/stop_stop/800_450'
inDir['NoPU']['1100_450'] = base+'/NoPU/ML/stop_stop/1100_450'
inDir['NoPU']['1125_25']  = base+'/NoPU/ML/stop_stop/1125_25'
# 50PU
inDir['50PU']['483_300']  = base+'/50PU/stop_stop/483_300'
inDir['50PU']['566_375']  = base+'/50PU/stop_stop/566_375'
inDir['50PU']['800_450']  = base+'/50PU/stop_stop/800_450'
inDir['50PU']['1100_450'] = base+'/50PU/stop_stop/1100_450'
inDir['50PU']['1125_25']  = base+'/50PU/stop_stop/1125_25'
# 140PU
inDir['140PU']['483_300']  = base+'/140PU/stop_stop/483_300'
inDir['140PU']['566_375']  = base+'/140PU/stop_stop/566_375'
inDir['140PU']['800_450']  = base+'/140PU/stop_stop/800_450'
inDir['140PU']['1100_450'] = base+'/140PU/stop_stop/1100_450'
inDir['140PU']['1125_25']  = base+'/140PU/stop_stop/1125_25'
#--------------------------------------------- end of sample properties

# 
from ROOT import TFile
from glob import glob
from sys import exit
def GetEntries(dirname):
	files = glob(dirname+'/*.root')
	if len(files)>1:
		print 'GetEntries: there is more than 1 root file in '+dirname
		exit(0)
	print dirname,files
	file=TFile(files[0])
	tree = file.Get("delphTree")
	return tree.GetEntries()
	
# do it
for scene in scenarios:
	for samp in samples:
		if do[scene][samp]: 
			f=''
			for i in range(len(dirsHT[samp])):
				entries = GetEntries(inDir[scene][samp]+dirsHT[samp][i])
				f=f+inDir[scene][samp]+dirsHT[samp][i]+' '+str(weights[samp][i]/entries)+' '
			print f,samp,scene
			if samp in noWeight:
				reader(f,scene+'_'+samp,False)
			else:   
				reader(f,scene+'_'+samp)





