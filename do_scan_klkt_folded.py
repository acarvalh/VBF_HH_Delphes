#!/usr/bin/env python
import os, sys, time,math
from optparse import OptionParser
parser = OptionParser()

parser.add_option("--model", type="string", dest="name", help="HHGG contact", default=0)
parser.add_option("--exp", type="string", dest="exp", help="HHGG contact", default=0)
parser.add_option("--bunch", type="int", dest="bunch", help="HHGG contact", default=0)
(options, args) = parser.parse_args()


#from HHStatAnalysis.AnalyticalModels.NonResonantModel import NonResonantModel
#import shutil,subprocess
#model = NonResonantModel()
#exp = "ATLAS"
exp = options.exp # "CMS"
name = options.name # "-klkt"
#name = "-2HDM"
#name = "-SingExpZ2"
#name = "-CompTrip"

if exp =="ATLAS" : nev = 3
if exp =="CMS" : nev = 6


if ("-CompTrip" in options.name) :
   # sinbeta
   par1min = 1
   par1max = 500
   par1step = 5
   # m = 21 bunches
   par2min = 200
   par2max = 3001
   par2step = 100
if  ("-RealTrip" in options.name): 
   # sinbeta
   par1min = 1
   par1max = 500
   par1step = 5
   # m 
   par2min = 200
   par2max = 3001
   par2step = 100
if "-2HDM" in options.name : 
   # sinbeta
   par1min = -100
   par1max = 101
   par1step = 10
   # Z6 
   par2min = -20
   par2max = 21
   par2step = 2
   mH = 300
if "-SingExpZ2" in options.name:
   # cosalpha
   par1min = 60
   par1max = 101
   par1step = 5
   # kl 
   par2min = -80
   par2max = 81
   par2step = 10
if "-SingSpoZ2" in options.name:
   # cosalpha
   par2min = 0
   par2max = 101
   par2step = 5
   # one time 
   par1min = 0
   par1max = 1
   par1step = 1
if "-VLQT" in options.name:
   # cosalpha
   par2min = 0
   par2max = 51
   par2step = 2
   # one time 
   par1min = 0
   par1max = 1
   par1step = 1
if "-VLLE" in options.name:
   # cosalpha
   par2min = 0
   par2max = 31
   par2step = 1
   # one time 
   par1min = 0
   par1max = 1
   par1step = 1
if "-c2kt" in options.name:
   # c2
   par2min = -70
   par2max = 70
   par2step = 5
   # one time 
   par1min = 50
   par1max = 26
   par1step = 2
if "-c2cg" in options.name:
   # c2
   par2min = -70
   par2max = 70
   par2step = 5
   # one time 
   par1min = 40
   par1max = 40
   par1step = 1


else : print "model not implemented"



v = 246.
lamSM = 0.012
ini=int(par1min+options.bunch*(par1step))

print (ini, par1min,par1max,par1step," ",par2min,par2max,par2step)

for  par1 in range(ini,ini+par1step,par1step) :
  print par1
  for par2 in range(par2min,par2max,par2step) : # 26
    if "-CompTrip" in options.name : 
      sinbeta = float(par1/1000.0) 
      m =  float(par2*1.0)  
      fpar1 = sinbeta
      fpar2 = m
      kt = 1 - 2*math.pow(sinbeta,2)
      kl = 1 + (2*math.pow(sinbeta,2)*(3 + (4*math.pow(m,2))/(math.pow(v,2)*lamSM)))
      c2 = -4*math.pow(sinbeta,2)
      cg = 0
      c2g = 0
    if "-RealTrip" in options.name : 
      sinbeta = float(par1/1000.0) 
      m =  float(par2*1.0)  
      fpar1 = sinbeta
      fpar2 = m
      kt = 1 
      kl = 1 + (4*math.pow(sinbeta,2)*(3 + math.pow(m,2)/(math.pow(v,2)*lamSM)))
      c2 = -2*math.pow(sinbeta,2)
      cg = 0
      c2g = 0
    if "-2HDM" in options.name : 
      cosbeta = float(par1/100.0) 
      z6 =  float(par2*1.0)
      fpar1 = cosbeta
      fpar2 = z6  
      tan = math.tan(math.acos(cosbeta))
      kt =  1 - math.pow(v,2)/math.pow(mH,2)*z6/tan
      kl = 1 - math.pow(v,2)/math.pow(mH,2)*1.5*z6^2/(2*lamSM)
      c2 = -(math.pow(v,2)/math.pow(mH,2))*1.5*z6/tan 
      cg = 0
      c2g = 0
    if "-SingExpZ2" in options.name:
      cosbeta = float(par1/100.0) 
      linear =  float(par2*1.0)  
      fpar1 = cosbeta
      fpar2 = linear
      tan2=math.pow(math.tan(math.acos(cosbeta)),2)
      kt = 1 - tan2/2
      kl = 1 - 1.5*tan2 + tan2*linear
      c2 = -tan2/2
      cg = 0
      c2g = 0
    if "-SingSpoZ2" in options.name:
      cosbeta = float(par2/100.0) 
      fpar2 = cosbeta
      fpar1 = 1.0
      tan2=math.pow(math.tan(math.acos(cosbeta)),2)
      kt = 1 - tan2/2
      kl = 1 - 1.5*tan2
      c2 = -tan2/2
      cg = 0
      c2g = 0
    if "-VLQT" in options.name:
      ratio = float(par2*0.0001) 
      fpar2 = ratio
      fpar1 = 1.0
      kt = 1 - pow(ratio*v,2)/2.0
      kl = 1.0
      c2 = -3*1 +pow(ratio*v,2)/4.0
      cg = 0
      c2g = 0
    if "-VLLE" in options.name:
      ratio = float(par2*0.001) 
      fpar2 = ratio
      fpar1 = 1.0
      kt = 1 + pow(ratio*v,2)/4.0
      kl = 1.0 + pow(ratio*v,2)/4.0
      c2 = 0
      cg = 0
      c2g = 0
    if "-c2kt" in options.name:
      fpar2 = float(par2/10.0)
      fpar1 = float(par1/10.0)
      kt = fpar2
      if "-c2kt_kl15" in options.name: kl = 15.0
      elif "-c2kt_kl15" in options.name: kl = -15.0
      else : kl = 1.0
      c2 = fpar1
      cg = 0
      c2g = 0
    if "-c2cg" in options.name:
      fpar2 = float(par2/100.0)
      fpar1 = float(par1/10.0)
      kt = 1.0
      if "-c2cg_kl15" in options.name: kl = 15.0
      elif "-c2cg_kl15" in options.name: kl = -15.0
      else : kl = 1.0
      c2 = fpar1
      cg = fpar2
      c2g = -1*fpar2

    else : print "model not implemented"
    print (fpar1,fpar2," ", kl,kt,c2,cg,c2g)
    #"""
    if 1>0  : #neventsTH >nev/0.25  and neventsTH <nev/0.01 : # neventsTH >nev/0.15 and neventsTH <nev/0.01 : 
      print "pottential to probe"
      os.system('python recastHH.py --kl '+str(kl)+' --kt '+str(kt)+' --c2 '+str(c2)+' --cg '+str(cg)+' --c2g '+str(c2g)+' --exp '+exp+' --out '+name+' --bunch '+str(options.bunch)+\
                ' --par1 '+str(fpar1)+' --par2 '+str(fpar2))
    else : 
      print "not/too constrained"
      ftree = open("efficiency_scan_"+exp+name+".txt", 'a+') 
      ftree.write(str(s))
      ftree.write(' ')
      ftree.write(str(kt))
      ftree.write(' ')
      ftree.write(str(c2))
      ftree.write(' ')
      ftree.write(str(0.0))
      ftree.write(' ')
      ftree.write(str(0.0))
      ftree.write(' ')
      if neventsTH >= nev/0.025 : ftree.write(str(float(10000.0)))
      else : ftree.write(str(float(0.0)))
      ftree.write('\n')
      ftree.close()
    #"""
