import numpy as np
from ROOT import TFile, gROOT, gStyle , TLatex
from rootpy.plotting import Hist, Hist2D
from rootpy.plotting import Canvas, Legend , Pad 
from rootpy.io import root_open
import csv
import json
import itertools as it
from collections import defaultdict
import scipy
from copy import deepcopy


#! /usr/bin/env python
# Analytical reweighting implementation for H->4b
# This file is part of https://github.com/cms-hh/HHStatAnalysis.
# python nonResonant_test_v0.py --kl 1 --kt 1 

# compiling
from optparse import OptionParser
import ROOT
import numpy as np
from HHStatAnalysis.AnalyticalModels.NonResonantModel import NonResonantModel


parser = OptionParser()
parser.add_option("--kl", type="float", dest="kll", help="Multiplicative factor in the H trilinear wrt to SM")
parser.add_option("--kt", type="float", dest="ktt", help="Multiplicative factor in the H top Yukawa wrt to SM")

parser.add_option("--c2", type="float", dest="c22", help="ttHH with triangle loop", default=0)
parser.add_option("--cg", type="float", dest="cgg", help="HGG contact", default=0)
parser.add_option("--c2g", type="float", dest="c2gg", help="HHGG contact", default=0)
parser.add_option("--doPlot", action='store_true', default=False, dest='doPlot', 
                    help='calculate the limit in the benchmark poin specified')
parser.add_option("--out", type="string", dest="out", help="HHGG contact", default=0)
parser.add_option("--exp", type="string", dest="exp", help="HHGG contact", default=0)
parser.add_option("--bunch", type="int", dest="bunch", help="HHGG contact", default=0)
parser.add_option("--par1", type="float", dest="par1", help="HHGG contact", default=-10000)
parser.add_option("--par2", type="float", dest="par2", help="HHGG contact", default=-10000)

(options, args) = parser.parse_args()
print " "
kl = options.kll
kt = options.ktt
c2 = options.c22
cg = options.cgg
c2g = options.c2gg

#print "events for V0 (the same of the fullsim version of Moriond 2016) \n We sum SM + box + the benchmarks from 2-13"
#if c2 != 0 or cg != 0 or c2g != 0 :  print "The analytical function is not yet implemented"

###########################################################
# read events and apply weight
###########################################################
def main(kl,kt,c2,cg,c2g):
  print (kl,kt,c2,cg,c2g)
  #if 1 > 0 :
  # declare the 2D ===> should be global variable
  model = NonResonantModel()
  # obtaining BSM/SM coeficients
  dumb = model.ReadCoefficients("../../HHStatAnalysis/AnalyticalModels/data/coefficientsByBin_extended_3M_costHHSim_19-4.txt") 
  counteventSM=0
  sumWeight=0
  #neventsTH = model.functionGF(kl,kt,c2,cg,c2g,model.A13tev)*3.2*33.45*0.0026
  #if neventsTH/0.4 <3 or neventsTH/0.001 >4 :
  #  print "not/too constrained" 
  #  return 
  # print "sum of weights calculated" , calcSumOfWeights 
  # read the events
  file=ROOT.TFile("/afs/cern.ch/work/a/acarvalh/HH4b_HL-LHC/VBF_HH_Delphes/"+str(options.exp)+"trees/"+str(options.exp)+"_HH_GF-SUM_JHEP.root")
  tree=file.Get("TCVARS6")
  nev = tree.GetEntries()
  # declare the histograms 
  CalcMhh = np.zeros((nev))
  CalcCost = np.zeros((nev))
  CalcPtHgg = np.zeros((nev))
  CalcPtHH = np.zeros((nev))
  CalcWeight = np.zeros((nev))
  CalcMhhReco = np.zeros((nev))
  CalcMXReco = np.zeros((nev))
  countevent = 0
  ####################
  histfilename="../Support/NonResonant/Distros_5p_SM3M_sumBenchJHEP_13TeV_19-4.root"
  #histtitle= "H1bin3" #"SumV0_AnalyticalBinExt"
  histtitle= "H1bin4" #"SumV0_AnalyticalBinExt"
  bm = model.getCluster(kl, kt,c2,cg,c2g,histfilename,histtitle)
  print "Closest benchmark is: ",  bm 
  fileHH=ROOT.TFile(histfilename)
  sumHAnalyticalBin = fileHH.Get(histtitle)
  calcSumOfWeights = model.getNormalization(kl, kt,c2,cg,c2g,histfilename,histtitle)  # this input is flexible, tatabb may have only the SM
  #for kll in range(-5,5) : model.getNormalization(kll, kt,sumHBenchBin)
  for iev in range(0,nev) :
      tree.GetEntry(iev)
      mhh = tree.MHH_Genlevel
      cost = tree.Costheta_Genlevel
      ptHgg = tree.pTgg_AllSel
      ptHH = tree.pThh_AllSel
      # to add cut
      mHHreco = tree.Mhh_AllSel
      mggreco = tree.Mgg_AllSel
      mXreco = tree.Mx_AllSel
      mhhcost= [mhh,cost,ptHgg,ptHH] # to store [mhh , cost] of that event
      # find the Nevents from the sum of events on that bin
      bmhh = sumHAnalyticalBin.GetXaxis().FindBin(mhh)
      bcost = sumHAnalyticalBin.GetYaxis().FindBin(abs(cost))
      effSumV0 = sumHAnalyticalBin.GetBinContent(bmhh,bcost)  # quantity of simulated events in that bin (without cuts)
      weight = model.getScaleFactor(mhh , cost,kl, kt,c2,cg,c2g, effSumV0 , calcSumOfWeights)
      #############################################
      # fill histograms to test
      #############################################
      if weight > 0: 
               #print countevent
               CalcMhh[countevent] = float(mhhcost[0]) 
               CalcCost[countevent] = float(abs(mhhcost[1])) 
               CalcPtHgg[countevent] = float(mhhcost[2]) 
               CalcPtHH[countevent] = float(mhhcost[3]) 
               CalcWeight[countevent] = weight 
               CalcMhhReco[countevent] = float(mHHreco) 
               CalcMXReco[countevent] = float(mXreco) 
               countevent+=1
               sumWeight+=weight
      else : " negative weight "
  print "plotted histogram reweighted from ",countevent," events, ", float(100*(nev-countevent)/nev)," % of the events was lost in empty bins in SM simulation"
  print "sum of weights",sumWeight," correction",calcSumOfWeights 
  ftree = open("efficiency_scan_"+str(options.exp)+str(options.out)+"_bunch_"+str(options.bunch)+".txt", 'a+') 
  if options.par1 != -10000 and options.par2 != -10000 :
    ftree.write(str(options.par1))
    ftree.write(' ')
    ftree.write(str(options.par2))
    ftree.write(' ')    
  ftree.write(str(kl))
  ftree.write(' ')
  ftree.write(str(kt))
  ftree.write(' ')
  ftree.write(str(c2))
  ftree.write(' ')
  ftree.write(str(cg))
  ftree.write(' ')
  ftree.write(str(c2g))
  ftree.write(' ')
  ftree.write(str(float(sumWeight)))
  ftree.write(' ')
  ftree.write(str(bm))
  ftree.write('\n')
  ############################################################################################################################
  # Draw test histos
  ###############################################
  drawtest =-1
  klJHEP=[1.0, 7.5,  1.0,  1.0,  -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0]
  ktJHEP=[1.0, 1.0,  1.0,  1.0,  1.5,  1.0, 1.0, 1.0, 1.0,  1.0, 1.5,  1.0, 1.0]
  c2JHEP=[0.0, -1.0, 0.5, -1.5, -3.0,  0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0]
  cgJHEP=[0.0, 0.0, 0.6,  0.0, 0.0,   0.8, 0.2, 0.2, -1.0, -0.6, 0.0, 1.0, 0.0]
  c2gJHEP=[0.0, 0.0, 1.0, -0.8, 0.0, -1.0, -0.2,-0.2,  1.0,  0.6, 0.0, -1.0, 0.0]
  # python recastHH.py --kl 7.5 --kt 1 --c2 -1
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 0.5 --cg -0.8 --c2g 0.6
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 -1.5 --cg 0.0 --c2g -0.8
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 0.0 --cg 0.8 --c2g -1.0
  # python recastHH.py --kl 1.0 --kt 1.0 --c2 1.0 --cg -0.6 --c2g 0.6
  for sam in range(0,13) :
    if kl == klJHEP[sam] and kt == ktJHEP[sam] and c2 ==c2JHEP[sam] and cg == cgJHEP[sam] and c2g ==c2gJHEP[sam] : 
       fileteste=ROOT.TFile("/afs/cern.ch/work/a/acarvalh/HH4b_HL-LHC/VBF_HH_Delphes/CMStrees/"+str(options.exp)+"_HH_GF-BM_"+str(sam)+".root")
       treeteste=fileteste.Get("TCVARS6")
       drawtest = sam
  ###########################################################################
  if drawtest>-1 : 
     nevtest = treeteste.GetEntries()
     CalcMhhTest = np.zeros((nevtest))
     CalcCostTest = np.zeros((nevtest))
     CalcPtHggTest = np.zeros((nevtest))
     CalcPtHHTest = np.zeros((nevtest))
     CalcMhhRecoTest = np.zeros((nevtest))
     CalcMXRecoTest = np.zeros((nevtest))
     countevent=0
     print drawtest,"sum of weights",sumWeight," correction",calcSumOfWeights," plain eff ", float(float(nevtest)/100000)
     for iev in range(0,nevtest) :
       treeteste.GetEntry(iev)
       ######
       mHHtest = treeteste.MHH_Genlevel
       costtest = treeteste.Costheta_Genlevel
       ptHggtest = treeteste.pTgg_AllSel
       ptHHtest = treeteste.pThh_AllSel
       # to add cut
       mHHrecotest = treeteste.Mhh_AllSel
       mggrecotest = treeteste.Mgg_AllSel
       mXrecotest = treeteste.Mx_AllSel
       mhhcosttest= [mHHtest,costtest,ptHggtest,ptHHtest] # to store [mhh , cost] of that event
       #
       CalcMhhTest[countevent] = float(mhhcosttest[0]) 
       CalcCostTest[countevent] = float(abs(mhhcosttest[1])) 
       CalcPtHggTest[countevent] = float(mhhcosttest[2]) 
       CalcPtHHTest[countevent] = float(mhhcosttest[3])
       CalcMhhRecoTest[countevent] = float(mHHrecotest)  
       CalcMXRecoTest[countevent] = float(mXrecotest)  
       countevent+=1
  else : 
     CalcMhhTest = np.zeros((1))
     CalcCostTest = np.zeros((1))
     CalcPtHggTest = np.zeros((1))
     CalcPtHHTest = np.zeros((1))
     CalcMhhRecoTest = np.zeros((1))
     CalcMXRecoTest = np.zeros((1))
  ############################################################################################################################ 
  if options.doPlot : 
    model.plotting(kl,kt,c2,cg,c2g,CalcMhh,CalcCost,CalcPtHgg,CalcPtHH,CalcMhhReco,CalcMXReco,CalcWeight,CalcMhhTest,CalcCostTest,CalcPtHggTest,CalcPtHHTest,CalcMhhRecoTest,CalcMXRecoTest,drawtest)

##########################################
if __name__ == "__main__": 
  #for  s in range(-40,-39) : #,40,2) :
  #  for t in range(-30,-29) :  #31) : 
  #    kl=s
  #    kt=float(t/5)
  #    c2=0.0
  #    cg=0.0
  #    c2g=0.0
  #    main(kl,kt,c2,cg,c2g)
  main(kl,kt,c2,cg,c2g)
