#!/usr/bin/env python
# to run: ./CMS_HH_analysis_code_GF-BM.py /eos/user/a/aaggarwa/delphes_GF-BM/CMS_delphes_GF_HH_1.root 1
import os, sys, time,math
import ROOT
#from ROOT import TLatex,TPad,TList,TH1,TH1F,TH2F,TH1D,TH2D,TFile,TTree,TCanvas,TLegend,TBranch,SetOwnership,gDirectory,TObject,gStyle,gROOT,TLorentzVector,TGraph,TMultiGraph,TColor,TAttMarker,TLine,TDatime,TGaxis,TF1,THStack,TAxis,TStyle,TPaveText,TAttFill,TF2, gPad, TGaxis, TChain,TClass 
from array import array
import numpy as np
pow = ROOT.TMath.Power
import bisect
from optparse import OptionParser
# Delphes headers
ROOT.gInterpreter.Declare('#include "ExRootAnalysis/ExRootTreeReader.h"')
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gSystem.Load("libDelphes")
# Fastjet headers = just to sort by PT
#ROOT.fastjet
#from ROOT.fastjet import *
#ROOT.gInterpreter.Declare('#include "fastjet/PseudoJet.hh"')
# need to fix the path there
#import ExRootAnalysis
#from ExRootAnalysis import ExRootTreeReader
parser = OptionParser()
if len(sys.argv) < 2:
    print " Usage: Example1.py input_file"
    sys.exit(1)
inputFile = sys.argv[1]
outputFile = sys.argv[2]
#############################################################
# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)
# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()
print "The tree have "+str(numberOfEntries)+" events"
# Get pointers to branches used in this analysis
branchEvent = treeReader.UseBranch("Event")
branchJet = treeReader.UseBranch("Jet")
branchParticle = treeReader.UseBranch("Particle")
branchPhoton = treeReader.UseBranch("Photon")
negative=0
#############################################################
# Declare histograms
#############################################################

Eff = 0.0 
Nevents = 0.0
############################################################# 
f1 = ROOT.TFile("CMStrees/CMS_HH_GF-BM_" + str(outputFile)+".root", 'RECREATE')
o_Mgg_PhSel = np.zeros(1, dtype=float)
o_pTgg_PhSel = np.zeros(1, dtype=float)
o_Mjj_JetSel = np.zeros(1, dtype=float)
o_pTjj_JetSel = np.zeros(1, dtype=float)
o_pTj1_JetSel = np.zeros(1, dtype=float)
o_pTj2_JetSel = np.zeros(1, dtype=float)
o_Eta_JetSel = np.zeros(1, dtype=float)
o_ExtraJets_JetSel = np.zeros(1, dtype=float)
o_Eta_noHjj_JetSel = np.zeros(1, dtype=float)
o_MaxDelEta_NoHjj_JetSel = np.zeros(1, dtype=float)
o_Mgg_AllSel = np.zeros(1, dtype=float)
o_pTgg_AllSel = np.zeros(1, dtype=float)
o_Mjj_AllSel = np.zeros(1, dtype=float)
o_pTjj_AllSel = np.zeros(1, dtype=float)
o_pTj1_AllSel = np.zeros(1, dtype=float)
o_pTj2_AllSel = np.zeros(1, dtype=float)
o_ExtraJets_AllSel = np.zeros(1, dtype=float)
o_Eta_AllSel = np.zeros(1, dtype=float)
o_Eta_noHjj_AllSel = np.zeros(1, dtype=float)
o_MaxDelEta_NoHjj_AllSel = np.zeros(1, dtype=float)
o_Mhh_AllSel = np.zeros(1, dtype=float)
o_costheta_AllSel = np.zeros(1, dtype=float)
o_pThh_AllSel = np.zeros(1, dtype=float)
o_Mx_AllSel = np.zeros(1, dtype=float)
o_Mhh_GenLevel = np.zeros(1, dtype=float)
o_costheta_GenLevel = np.zeros(1, dtype=float)

o_Mhh_GenLevelnocut = np.zeros(1, dtype=float)
o_costheta_GenLevelnocut = np.zeros(1, dtype=float)

t = ROOT.TTree("TCVARS", "Limit tree for HH->bbgg analyses")
t0 = ROOT.TTree("TCVARS0", "nocuts")
t1 = ROOT.TTree("TCVARS1", "Limit tree for HH->bbgg analyses")
t2 = ROOT.TTree("TCVARS2", "Limit tree for HH->bbgg analyses")
t3 = ROOT.TTree("TCVARS3", "Limit tree for HH->bbgg analyses")
t4 = ROOT.TTree("TCVARS4", "Limit tree for HH->bbgg analyses")
t5 = ROOT.TTree("TCVARS5", "Limit tree for HH->bbgg analyses")
t6 = ROOT.TTree("TCVARS6", "Limit tree for HH->bbgg analyses")
t7 = ROOT.TTree("TCVARS7", "Limit tree for HH->bbgg analyses")
t8 = ROOT.TTree("TCVARS8", "Limit tree for HH->bbgg analyses")
#t9 = ROOT.TTree("TCVARS9", "Limit tree for HH->bbgg analyses")
#t10 = ROOT.TTree("TCVARS10", "Limit tree for HH->bbgg analyses")
t.Branch("Mgg_PhSel", o_Mgg_PhSel, "o_Mgg_PhSel/D");
t.Branch("pTgg_PhSel", o_pTgg_PhSel, "o_pTgg_PhSel/D");
t0.Branch("MHH_Genlevelnocut", o_Mhh_GenLevelnocut, "o_Mhh_GenLevelnocut/D")
t0.Branch("Costheta_Genlevelnocut", o_costheta_GenLevelnocut, "o_costheta_GenLevelnocut/D")
t1.Branch("Mjj_JetSel", o_Mjj_JetSel, "o_Mjj_JetSel/D");
t1.Branch("pTjj_JetSel", o_pTjj_JetSel, "o_pTjj_JetSel/D");
t1.Branch("pTj1_JetSel", o_pTj1_JetSel, "o_pTj1_JetSel/D");
t1.Branch("pTj2_JetSel", o_pTj2_JetSel, "o_pTj2_JetSel/D");
t2.Branch("Eta_JetSel", o_Eta_JetSel, "o_Eta_JetSel/D");
t3.Branch("ExtraJets_JetSel", o_ExtraJets_JetSel, "o_ExtraJets_JetSel/D");
t5.Branch("Eta_noHjj_JetSel", o_Eta_noHjj_JetSel, "o_Eta_noHjj_JetSel/D");
t4.Branch("MaxDelEta_NoHjj_JetSel", o_MaxDelEta_NoHjj_JetSel, "o_MaxDelEta_NoHjj_JetSel/D");
t6.Branch("Mgg_AllSel", o_Mgg_AllSel, "o_Mgg_AllSel/D");
t6.Branch("pTgg_AllSel", o_pTgg_AllSel, "o_pTgg_AllSel/D");
t6.Branch("Mjj_AllSel", o_Mjj_AllSel, "o_Mjj_AllSel/D");
t6.Branch("pTjj_AllSel", o_pTjj_AllSel, "o_pTjj_AllSel/D");
t6.Branch("pTj1_AllSel", o_pTj1_AllSel, "o_pTj1_AllSel/D");
t6.Branch("pTj2_AllSel", o_pTj2_AllSel, "o_pTj2_AllSel/D");
t7.Branch("Eta_AllSel", o_Eta_AllSel, "o_Eta_AllSel/D");
t6.Branch("ExtraJets_AllSel", o_ExtraJets_AllSel, "o_ExtraJets_AllSel/D");
t8.Branch("Eta_noHjj_AllSel", o_Eta_noHjj_AllSel, "o_Eta_noHjj_AllSel/D");
t6.Branch("MaxDelEta_NoHjj_AllSel", o_MaxDelEta_NoHjj_AllSel, "o_MaxDelEta_NoHjj_AllSel/D");
t6.Branch("Mhh_AllSel", o_Mhh_AllSel, "o_Mhh_AllSel/D");
t6.Branch("Costheta_AllSel", o_costheta_AllSel, "o_Costheta_AllSel/D"); # to calculate
t6.Branch("pThh_AllSel", o_pThh_AllSel, "o_pThh_AllSel/D");
t6.Branch("Mx_AllSel", o_Mx_AllSel, "o_Mx_AllSel/D");
t6.Branch("MHH_Genlevel", o_Mhh_GenLevel, "o_Mhh_GenLevel/D")
t6.Branch("Costheta_Genlevel", o_costheta_GenLevel, "o_costheta_GenLevel/D")

ftree = open("CMStrees/CMS_HH_GF-BM_" + str(outputFile)+"_events.txt", 'a+') 
ftree.write('o_Mgg_AllSel')
ftree.write(' ')
ftree.write('o_pTgg_AllSel')
ftree.write(' ')
ftree.write('o_Mjj_AllSel')
ftree.write(' ')
ftree.write('o_pTjj_AllSel')
ftree.write(' ')
ftree.write('o_pTj1_AllSel')
ftree.write(' ')
ftree.write('o_pTj2_AllSel')
ftree.write(' ')
ftree.write('o_ExtraJets_AllSel')
ftree.write(' ')
ftree.write('o_Eta_AllSel')
ftree.write(' ')
ftree.write('o_Eta_noHjj_AllSel')
ftree.write(' ')
ftree.write('o_MaxDelEta_NoHjj_AllSel')
ftree.write(' ')
ftree.write('o_Mhh_AllSel')
ftree.write(' ')
ftree.write('o_pThh_AllSel')
ftree.write(' ')
ftree.write('o_Mx_AllSel')
ftree.write(' ')
ftree.write('o_Mhh_GenLevel')
ftree.write(' ')
ftree.write('o_costheta_GenLevel')
ftree.write('\n')
counter = 0
counter2=0
for entry in range(0, treeReader.GetEntries()):
    #for entry in range(0, 500):
    counter+=1
    if counter ==1000 : 
       counter = 0
       counter2+=1
       print ('processed ',counter2,'k events')
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)
    weight = branchEvent.At(0).Weight
    if(weight < 0): 
        negative+=1
        continue
    
    ########################
    # GEN LEVEL PARTICLES
    ########################
    Higgses = []
    SUM =[]

    for part in range(0, branchParticle.GetEntries()):

       genparticle =  branchParticle.At(part)
       dumb100 = ROOT.TLorentzVector()
       dumb100.SetPxPyPzE(genparticle.Px,genparticle.Py,genparticle.Pz,genparticle.E)
       pdgCode = genparticle.PID
       IsPU = genparticle.IsPU
       status = genparticle.Status
       Pt = genparticle.PT
       Px = genparticle.Px
       Py = genparticle.Py
       Pz = genparticle.Pz
       Energy = genparticle.E
       Mass = genparticle.Mass

       if(IsPU == 0 and (pdgCode == 25)): # and status==statusH and abs(motherPID)  > 600000): 
          mother =  branchParticle.At(genparticle.M1)
          motherPID = mother.PID
          if(branchParticle.At(genparticle.D1).PID != 25 and (branchParticle.At(genparticle.D1).Status == 23 or branchParticle.At(genparticle.D1).Status == 1) and (branchParticle.At(genparticle.D2).Status == 23 or branchParticle.At(genparticle.D2).Status == 1)) :
	     Higgses.append(dumb100)          

    if (len(Higgses) >1):
                SUM = ROOT.TLorentzVector()
                SUM.SetPxPyPzE(Higgses[0].Px() + Higgses[1].Px(), Higgses[0].Py() + Higgses[1].Py(), Higgses[0].Pz() + Higgses[1].Pz() , Higgses[0].E() + Higgses[1].E())
                Mhh_GenLevel =SUM.M()

                P1boost = Higgses[0]
                P1boost.Boost(-SUM.BoostVector())
                costheta_GenLevel =  float(P1boost.CosTheta())
		o_Mhh_GenLevelnocut[0] = Mhh_GenLevel
		o_costheta_GenLevelnocut[0] = costheta_GenLevel

		o_Mhh_GenLevel[0] = Mhh_GenLevel
		o_costheta_GenLevel[0] = costheta_GenLevel

                t0.Fill()
    #else : print "not 2 higgses"
    #####################
    # PHOTONS
    #####################
    RecoPhotons_NoCut =[]
    RecoPhotons =[]
    RecoPhotons_afterSel =[]
    TwoPhotons_afterAllSel = False
    Mgg = 0.0
    pTgg = 0.0
    for part in range(0, branchPhoton.GetEntries()):
        photon =  branchPhoton.At(part)
        dumb1 = ROOT.TLorentzVector()
        dumb1.SetPtEtaPhiE(photon.PT,photon.Eta,photon.Phi,photon.E)
	RecoPhotons_NoCut.append(dumb1)
	if (abs (photon.Eta) <  2.5):
        	RecoPhotons.append(dumb1)
    if (len(RecoPhotons)> 1):
    	if ( RecoPhotons[0].Pt() > 30 and RecoPhotons[1].Pt() > 20 and ( (RecoPhotons[0].Pt())/((RecoPhotons[0] + RecoPhotons[1]).M()) ) > 1/3 and ( (RecoPhotons[1].Pt())/((RecoPhotons[0] + RecoPhotons[1]).M()) ) > 1/4 and (RecoPhotons[0] + RecoPhotons[1]).M() > 100 and (RecoPhotons[0] + RecoPhotons[1]).M() < 180 ):
		RecoPhotons_afterSel = RecoPhotons
		TwoPhotons_afterAllSel = True
		Mgg = (RecoPhotons[0] + RecoPhotons[1]).M()
		pTgg = (RecoPhotons[0] + RecoPhotons[1]).Pt()
		o_Mgg_PhSel[0] = Mgg
		o_pTgg_PhSel[0] = pTgg
		t.Fill()
    ##############################################################################################################
    #####################
    # Jets
    #####################
    RecoJets =[]
    RecoJets_after_FirstSel =[]
    RecoJets_MedBtag =[]
    RecoVBFJets =[]
    RecoVBFJets_noHjj =[]
    TwoJets_afterAllSel = False  
 
    for part1 in range(0, branchJet.GetEntries()):
        jet =  branchJet.At(part1)
        dumb2 = ROOT.TLorentzVector()
        dumb2.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
        RecoJets.append(dumb2)
	count = 0
    	for part2 in range(0, len(RecoPhotons_NoCut)):
		if (RecoPhotons_NoCut[part2].DeltaR(dumb2) < 0.4):
			count = count +1
	if (jet.PT > 25 and abs(jet.Eta) < 5 and count == 0 ):
		RecoVBFJets.append(dumb2)
	if (jet.PT > 25 and abs(jet.Eta) < 2.4 and count == 0):
		RecoJets_after_FirstSel.append(dumb2)
     		if (jet.BTag & 1 << 0):
			RecoJets_MedBtag.append(dumb2)
    #print (len(RecoPhotons),len(RecoJets_after_FirstSel),len(RecoJets_MedBtag))
    A1 =0
    A2 =1
    Mjj = 0.0
    pTjj = 0.0
    pTj1 = 0.0
    pTj2 = 0.0
   
    if (len(RecoJets_MedBtag)> 1):
	    for q in range(0, len(RecoJets_MedBtag)-1):
		for r in range(q+1, len(RecoJets_MedBtag)):
			if(((RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2]).Pt()) < ((RecoJets_MedBtag[q] + RecoJets_MedBtag[r]).Pt())):
				A1 = q
				A2 = r
	    if((RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2]).M() < 200 and (RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2]).M() > 80):
			TwoJets_afterAllSel= True
			Mjj = (RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2]).M()
                        pTjj = (RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2]).Pt()
			pTj1 = (RecoJets_MedBtag[A1]).Pt()
			pTj2 = (RecoJets_MedBtag[A2]).Pt()
			o_Mjj_JetSel[0] = Mjj
			o_pTjj_JetSel[0] = pTjj
			o_pTj1_JetSel[0] = pTj1
			o_pTj2_JetSel[0] = pTj2
			t1.Fill()

    if(TwoJets_afterAllSel == True):
	    for i1 in range(0, len(RecoVBFJets)):
		o_Eta_JetSel[0] = RecoVBFJets[i1].Eta()
		t2.Fill()

	    	if (RecoVBFJets[i1] != RecoJets_MedBtag[A1] and RecoVBFJets[i1] != RecoJets_MedBtag[A2]):
			dumb3 = ROOT.TLorentzVector()
			dumb3.SetPtEtaPhiM(RecoVBFJets[i1].Pt(),RecoVBFJets[i1].Eta(),RecoVBFJets[i1].Phi(),RecoVBFJets[i1].M())
			RecoVBFJets_noHjj.append(dumb3)
	    o_ExtraJets_JetSel[0] = len(RecoVBFJets_noHjj)
	    t3.Fill()
       
    #print(TwoJets_afterAllSel,TwoPhotons_afterAllSel)
    ####### VBF jets selection for 2B category ########
    VBFjet_maxEtaPair_noHjj_1 =0
    VBFjet_maxEtaPair_noHjj_2 =1
    MaxDelEta =0.0

    if (len (RecoVBFJets_noHjj) >1):
	    for z111 in range(0, len(RecoVBFJets_noHjj)-1): 
		for z211 in range (z111+ 1,len(RecoVBFJets_noHjj) ) :  
			if (  (abs(RecoVBFJets_noHjj[VBFjet_maxEtaPair_noHjj_1].Eta() - RecoVBFJets_noHjj[VBFjet_maxEtaPair_noHjj_2].Eta())) < (abs(RecoVBFJets_noHjj[z111].Eta() - RecoVBFJets_noHjj[z211].Eta())) ):
				VBFjet_maxEtaPair_noHjj_1 = z111
				VBFjet_maxEtaPair_noHjj_2 = z211
	    MaxDelEta = (abs(RecoVBFJets_noHjj[VBFjet_maxEtaPair_noHjj_1].Eta() - RecoVBFJets_noHjj[VBFjet_maxEtaPair_noHjj_2].Eta()))
            o_MaxDelEta_NoHjj_JetSel[0] = MaxDelEta
	    t4.Fill()

    if (len (RecoVBFJets_noHjj) >0):
	    for u1 in range(0, len(RecoVBFJets_noHjj)): 
	        o_Eta_noHjj_JetSel[0] =  RecoVBFJets_noHjj[u1].Eta()
		t5.Fill()

    Mhh = 0.0
    pThh = 0.0
    Mx = 0.0

      ##################################################################################################   


    if(TwoPhotons_afterAllSel == True and TwoJets_afterAllSel == True):
	Mhh= (RecoPhotons[0] + RecoPhotons[1] + RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2] ).M()
	pThh = (RecoPhotons[0] + RecoPhotons[1] + RecoJets_MedBtag[A1] + RecoJets_MedBtag[A2] ).Pt()
	Mx = Mhh - Mgg - Mjj + 250
	o_Mgg_AllSel[0] = Mgg
	o_pTgg_AllSel[0] = pTgg
	o_Mjj_AllSel[0] = Mjj
	o_pTjj_AllSel[0] = pTjj
	o_pTj1_AllSel[0] = pTj1
	o_pTj2_AllSel[0] = pTj2
	o_ExtraJets_AllSel[0] = len(RecoVBFJets_noHjj)
	o_Mhh_AllSel[0] = Mhh
	o_pThh_AllSel[0] = pThh
	o_Mx_AllSel[0] = Mx
	#if (len(Higgses) >1):
	#	o_Mhh_GenLevel[0] = Mhh_GenLevel
	#	o_costheta_GenLevel[0] = costheta_GenLevel

        if (len (RecoVBFJets_noHjj) >1):
                o_MaxDelEta_NoHjj_AllSel[0] = MaxDelEta

	

        for d1 in range(0, len(RecoVBFJets)):
                o_Eta_AllSel[0] = RecoVBFJets[d1].Eta()
		t7.Fill()
	if (len (RecoVBFJets_noHjj) >0):
		for g1 in range(0, len(RecoVBFJets_noHjj)):
			o_Eta_noHjj_AllSel[0] = RecoVBFJets_noHjj[g1].Eta()
			t8.Fill()

	if (Mhh > 350):
		Nevents = Nevents +1 
                ftree.write(str(o_Mgg_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_pTgg_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_Mjj_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_pTjj_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_pTj1_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_pTj2_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_ExtraJets_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_Eta_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_Eta_noHjj_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_MaxDelEta_NoHjj_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_Mhh_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_pThh_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_Mx_AllSel[0]))
                ftree.write(' ')
                ftree.write(str(o_Mhh_GenLevel[0]))
                ftree.write(' ')
                ftree.write(str(o_costheta_GenLevel[0]))
                ftree.write('\n')
                t6.Fill()

ftree.close()        		
####################################################################

Eff = ((float)(Nevents) / (treeReader.GetEntries())) # 

print (str(outputFile),Eff)
print "saved efficiencies to the file: CMStrees/CMS_HH_GF-BM_" + str(outputFile)+"_eff.txt"
f = open("CMStrees/CMS_HH_GF-BM_" + str(outputFile)+"_eff.txt", 'a+') 
f.write(str(Eff))
f.write('\n')
f.close()
###########################################################################################
f1.Write()
f1.Close()
#############################################################################################
#############################################################################################			
