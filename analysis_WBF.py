#!/usr/bin/env python
# to run: ./analysis.py delphes_output.root
import os, sys, time,math
import ROOT
#from ROOT import TLatex,TPad,TList,TH1,TH1F,TH2F,TH1D,TH2D,TFile,TTree,TCanvas,TLegend,SetOwnership,gDirectory,TObject,gStyle,gROOT,TLorentzVector,TGraph,TMultiGraph,TColor,TAttMarker,TLine,TDatime,TGaxis,TF1,THStack,TAxis,TStyle,TPaveText,TAttFill,TF2, gPad, TGaxis, TChain,TClass 
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

#########################
# Cuts
#########################
Tau21cut = 1
Tau31cut = 1
PrunMass2 = 100
PrunMass3 = 1000
VLQresolution = 100
Hresolution=50
#########################
# categorization
#########################
VLQTag2 = 0
VLQTag1 = 0
HTag2 = 0
HTag1 = 0
HTag0 = 0
#
VLQTag2truth = 0
VLQTag1truth = 0
HTag2truth = 0
HTag1truth = 0
HTag0truth = 0
#
HRes1truth = 0
HRes0truth = 0
#########################
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
branchFatJet = treeReader.UseBranch("Jet")
branchParticle = treeReader.UseBranch("Particle")
branchPhoton = treeReader.UseBranch("Photon")
#############################################################
# Declare histograms
#############################################################
# genvariables
histNoOfPhoton = ROOT.TH1F("No of photons", "No of photons", 10, 0, 10)

histH1PT = ROOT.TH1F("H1_pt", "Leading H Gen P_{T}", 100, 0.0, 3000.0)
histH1Mass = ROOT.TH1F("H1_mass", "Leading Gen M_{H}", 100, 40.0, 140.0)
histH2PT = ROOT.TH1F("H2_pt", "Sub-leading H Gen P_{T}", 100, 0.0, 3000.0)
histH2Mass = ROOT.TH1F("H2_mass", "Sub-leading  Gen M_{H}", 100, 40.0, 140.0)

histT1PT = ROOT.TH1F("T1_pt", "Leading Q Gen P_{T}", 100, 0.0, 3000.0)
histT1Mass = ROOT.TH1F("T1_mass", "Leading Gen M_{Q}", 100, 400.0, 3000.0)
histT2PT = ROOT.TH1F("T2_pt", "Sub-leading Q Gen P_{T}", 100, 0.0, 3000.0)
histT2Mass = ROOT.TH1F("T2_mass", "Sub-leading Gen M_{Q}", 100, 400.0, 3000.0)

histGenJ1PT = ROOT.TH1F("GenJ1_pt", "Leading GenJ P_{T}", 100, 0.0, 3000.0)
histGenJ2PT = ROOT.TH1F("GenJ2_pt", "Sub-leading GenJ P_{T}", 100, 0.0, 3000.0)

histNJets = ROOT.TH1F("NJets", "# Jets", 21, -0.5, 20.5)
histNFatJets = ROOT.TH1F("NFatJets", "# Fat Jets", 11, -0.5, 10.5)
histNBJets = ROOT.TH1F("NBJets", "#B Jets", 21, -0.5, 20.5)

histJ1Mass = ROOT.TH1F("J1_mass", "Leading J Mass", 100, 0.0, 3000.0)
histFatJ1Mass = ROOT.TH1F("FatJ1_mass", "Leading FatJ Mass", 100, 0.0, 3000.0)
histFatJ1Tau3 = ROOT.TH1F("FatJ1_tau3", "Leading FatJ #tau_{31}", 10, 0.0, 1)
histFatJ1Tau2 = ROOT.TH1F("FatJ1_tau2", "Leading FatJ #tau_{21}", 10, 0.0, 1)

histH1RecoMass = ROOT.TH1F("H1_mass_Reco", "Leading M_{H}", 100, 40.0, 1400.0)
histH2RecoMass = ROOT.TH1F("H2_mass_Reco", "Sub-leading  M_{H}", 100, 40.0, 1400.0)
#############################################################
# Loop over all events
#############################################################
negative=0

for entry in range(0, numberOfEntries):
    no_of_photons =0
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)
    weight = branchEvent.At(0).Weight
    #print weight
    if(weight < 0 ) : 
        negative+=1
        continue
    #####################
    # Gen-level particles
    #####################
    Higgses = []
    Topone = []
    GenJets = []
    QQ = True
    if (QQ) :
      statusH = 52# 22#52
      statusT = 62 #62
      statusGenJ = 62
    #print branchParticle.GetEntries()
    for part in range(0, branchPhoton.GetEntries()):
       no_of_photons = no_of_photons +1
    print "no of phtons:", no_of_photons
    histNoOfPhoton.Fill(no_of_photons)
    for part in range(0, branchParticle.GetEntries()):
       genparticle =  branchParticle.At(part)
       pdgCode = genparticle.PID
       IsPU = genparticle.IsPU
       status = genparticle.Status 
       # check if it is the correct status (for QQ the last 25 is 52 and the last topone 62)
       #print " pdgid "+ str(pdgCode)+" status "+str(status)
       if(IsPU == 0 and (pdgCode == 25) and status==statusH ): 
          Higgses.append(genparticle)
          #print " pdgid "+ str(pdgCode)+" status "+str(status)
       if (IsPU == 0 and (abs(pdgCode) > 600000) and status==statusT ): 
          Topone.append(genparticle)
          #print " pdgid "+ str(pdgCode)+" status "+str(status)
    
    if (len(Higgses) ==2) :
       # sort by pt = weird!
       if (Higgses[1].PT > Higgses[0].PT): 
          genparticle = Higgses[0]
          Higgses[0] = Higgses[1]
          Higgses[1] = genparticle
       histH1PT.Fill(Higgses[0].PT )
       histH2PT.Fill(Higgses[1].PT )
       histH1Mass.Fill(Higgses[0].Mass )
       histH2Mass.Fill(Higgses[1].Mass )
    else : print "not two Higggses " + str(len(Higgses))
    if (len(Topone) ==2) : # if there are 2
       if (Topone[1].PT > Topone[0].PT): 
          genparticle = Topone[0]
          Topone[0] = Topone[1]
          Topone[1] = genparticle    
       histT1PT.Fill(Topone[0].PT )
       histT2PT.Fill(Topone[1].PT )
       histT1Mass.Fill(Topone[0].Mass )
       histT2Mass.Fill(Topone[1].Mass )
    elif (len(Topone) ==1) :   
       histT1PT.Fill(Topone[0].PT )
       histT2PT.Fill(-100 )
       histT1Mass.Fill(Topone[0].Mass )
       histT2Mass.Fill(-100 )
    elif (len(Topone)>2 or len(Topone)<1) : print "too much (few) Quark partners"
    #####################
    # Jets - those are sorted by PT
    #####################
    #print "Jet properties"
    histNJets.Fill(branchJet.GetEntries())
    RecoJets = []
    RecoJetsBTag = []
    for part in range(0, branchJet.GetEntries()):
        jet =  branchJet.At(part)
        dumb = ROOT.TLorentzVector()
        dumb.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
        RecoJets.append(dumb)
        if (jet.BTag > 0 ) : RecoJetsBTag.append(dumb)
        #print "btag "+str(jet.BTag)
        #print "btag algo "+str(jet.BTagAlgo)
        if (part==0) : 
           histJ1Mass.Fill(jet.Mass)
    histNBJets.Fill(len(RecoJetsBTag))
    #####################
    # Fat Jets - those are sorted by PT
    #####################
    #print "FatJet properties"
    RecoFatJets3 = []
    RecoFatJets2 = []
    for part in range(0, branchFatJet.GetEntries()):
        jet =  branchFatJet.At(part)
        dumb = ROOT.TLorentzVector()
        dumb.SetPtEtaPhiM(jet.PT,jet.Eta,jet.Phi,jet.Mass)
        JMass = jet.Mass # it should be SoftDroppedP4
        if ( jet.Tau[3] > Tau31cut and JMass > PrunMass3 ) :  RecoFatJets3.append(dumb)
        if ( jet.Tau[2] > Tau21cut and JMass > PrunMass2 ) : RecoFatJets2.append(dumb)
        if (part==0) : 
            histFatJ1Mass.Fill(JMass) 
            histFatJ1Tau2.Fill(0)
            histFatJ1Tau3.Fill(0)
    histNFatJets.Fill(len(RecoFatJets2))
    ########################################################
    # Algoritm for categorization - add the mass selections
    ########################################################
    tolerance = 0.8
    if(len(RecoFatJets3)>1) :
       VLQTag2+=1
       # Doest it match Truth? - up to one
       count =0
       for ii in range(0,2): 
           dumb = ROOT.TLorentzVector()
           dumb.SetPtEtaPhiM(Topone[ii].Px,Topone[ii].Eta,Topone[ii].Phi,Topone[ii].Mass)
           for jj in range(0,2):
               #print str(Topone[ii].Mass)+" "+str(dumb.M()) +" "+ str(RecoFatJets3[jj].M())
               if ( RecoFatJets3[jj].DeltaR( dumb ) < tolerance) : count+=1
       if(count > 0) : VLQTag2truth+=1
       if(count > 2) : print "to much truth in 2 VLQ tag" 
       # fill histos
    elif(len(RecoFatJets3)>0 and len(RecoJetsBTag) > 1) :
       VLQTag1+=1
       # Doest it match Truth?
       count =0
       for ii in range(0,2): 
           dumb = ROOT.TLorentzVector()
           dumb.SetPtEtaPhiM(Topone[ii].Px,Topone[ii].Eta,Topone[ii].Phi,Topone[ii].Mass)
           if ( RecoFatJets3[0].DeltaR( dumb ) < tolerance) : count+=1
       if(count > 0) : VLQTag1truth+=1
       if(count > 1) : print "to much truth in 1 VLQ tag" 
       # fill histos
    elif(len(RecoFatJets2)>1) :
       HTag2+=1
       histH1RecoMass.Fill(RecoFatJets2[0].M())
       histH2RecoMass.Fill(RecoFatJets2[1].M())
       # Doest it match Truth? - only one is enought 
       count =0
       for ii in range(0,len(Higgses)): 
          dumb = ROOT.TLorentzVector()
          dumb.SetPtEtaPhiM(Higgses[ii].Px,Higgses[ii].Eta,Higgses[ii].Phi,Higgses[ii].Mass)
          for jj in range(0,2):
              if ( RecoFatJets2[jj].DeltaR( dumb ) < tolerance) : count+=1
       if(count > 0) : HTag2truth+=1 # print "to much truth, 2H "+str(count) 
       if(count > 2) : print "to much truth in 2 H tag" 
       # fill histos
    elif(len(RecoFatJets2)>0 and len(RecoJetsBTag) > 1) : # 
       HTag1+=1
       histH1RecoMass.Fill(RecoFatJets2[0].M())
       ######### Doest it match Truth?
       count =0
       for ii in range(0,len(Higgses)): 
           dumb = ROOT.TLorentzVector()
           dumb.SetPtEtaPhiM(Higgses[ii].Px,Higgses[ii].Eta,Higgses[ii].Phi,Higgses[ii].Mass)
           if ( RecoFatJets2[0].DeltaR( dumb ) < tolerance) : count+=1
       if(count > 0) : HTag1truth+=1 # print "to much truth, 1H "+str(count) 
       if(count > 1) : print "to much truth in 1 H tag " 
       ########## The other Higgs
       tomin = np.zeros((1000))
       tomini = []
       tominj = []
       count=0
       for ii in range(0,len(RecoJetsBTag)): 
           for jj in range(ii,len(RecoJetsBTag)): 
               #print str(count)+" "+str((RecoJetsBTag[ii]+RecoJetsBTag[jj]).M())+" "+str(RecoFatJets2[0].M())+" "+str(abs((RecoJetsBTag[ii]+RecoJetsBTag[jj]).M() - RecoFatJets2[0].M()))
               tomin[count] = abs((RecoJetsBTag[ii]+RecoJetsBTag[jj]).M() - RecoFatJets2[0].M())
               tomini.append(ii)
               tominj.append(jj)
               count+=1
       #print str(int(np.amin(tomin)))+" "+str((RecoJetsBTag[tomini[int(np.amin(tomin))]]+RecoJetsBTag[tominj[int(np.amin(tomin))]]).M())
       ResovedH = RecoJetsBTag[tomini[int(np.amin(tomin))]]+RecoJetsBTag[tominj[int(np.amin(tomin))]]
       histH2RecoMass.Fill(ResovedH.M())
       ########## The VLQs
       # fill histos
    else :
       HTag0+=1
       # Doest it match Truth?
       # fill histos
print "VLQTag2 VLQTag1 HTag2 HTag1 HTag0 Neg"
print str(VLQTag2)+" "+str(VLQTag1)+" "+str(HTag2)+" "+str(HTag1)+" "+str(HTag0)+" "+str(negative)
print str(VLQTag2truth)+" "+str(VLQTag1truth)+" "+str(HTag2truth)+" "+str(HTag1truth)+" "+str(HTag0truth)+" "+str(negative)

#######################
# save histos as image - test
#####################
c1=ROOT.TCanvas("c2","c2",200,50,600,600)
histH1PT.Draw()
c1.Print("histH1PT.png")
c1.Clear()

#######################
# save histos as root file - to superimpose
#####################
f = ROOT.TFile("histos_new.root", 'RECREATE')
histNoOfPhoton.Write()
histH1PT.Write()
histH2PT.Write()
histH1Mass.Write()
histH2Mass.Write()
histT1PT.Write()
histT2PT.Write()
histT1Mass.Write()
histT2Mass.Write()
histNJets.Write()
histNBJets.Write()
histNFatJets.Write()
histJ1Mass.Write()
histFatJ1Mass.Write()
histFatJ1Tau2.Write()
histFatJ1Tau3.Write()
histH1RecoMass.Write()
histH2RecoMass.Write()
#histGenJ1PT.Write()
#histGenJ2PT.Write()
f.Write()
f.Close()
