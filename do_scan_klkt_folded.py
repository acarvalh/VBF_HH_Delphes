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


if ("-CompTrip" in options.name) or ("-RealTrip" in options.name): 
   # sinbeta
   par1min = 0
   par1max = 100
   par1step = 5
   # m 
   par2min = 200
   par2max = 1000
   par2step = 100
if "-2HDM" in options.name : 
   # sinbeta
   par1min = -100
   par1max = 100
   par1step = 10
   # Z6 
   par2min = -20
   par2max = 20
   par2step = 2
   mH = 300
if "-SingExpZ2" :
   # cosalpha
   par1min = 0
   par1max = 100
   par1step = 5
   # kl 
   par2min = -50
   par2max = 50
   par2step = 2

else : print "model not implemented"

v = 246.
lamSM = 0.012
ini=int(par1min+options.bunch*(par1step+1))
for  par1 in range(par1min,par1max,par1step) :
  for par2 in range(par2min,par2max,par2step) : # 26
    if "-CompTrip" in options.name : 
      sinbeta = par1/100 
      m =  par2  
      kt = 1 - 2*math.pow(sinbeta,2)
      kl = 1 + (2*math.pow(sinbeta,2)*(3 + (4*math.pow(m,2))/(math.pow(v,2)*lamSM)))
      c2 = -4*math.pow(sinbeta,2)
      cg = 0
      c2g = 0
    if "-RealTrip" in options.name : 
      sinbeta = par1/100 
      m =  par2  
      kt = 1 
      kl = 1 + (4*math.pow(sinbeta,2)*(3 + math.pow(m,2)/(math.pow(v,2)*lamSM)))
      c2 = -2*math.pow(sinbeta,2)
      cg = 0
      c2g = 0
    if "-2HDM" in options.name : 
      cosbeta = par1/100 
      z6 =  par2  
      tan = math.tan(math.acos(cosbeta))
      kt =  1 - math.pow(v,2)/math.pow(mH,2)*z6/tan
      kl = 1 - math.pow(v,2)/math.pow(mH,2)*1.5*z6^2/(2*lamSM)
      c2 = -(math.pow(v,2)/math.pow(mH,2))*1.5*z6/tan 
      cg = 0
      c2g = 0
    if "-SingExpZ2" :
      cosbeta = par1/100 
      linear =  par2  
      tan2=math.pow(math.tan(math.acos(cosbeta)),2)
      kt = 1 - tan2/2
      kl = 1 - 1.5*tan2 + tan2*linear
      c2 = -tan2/2
      cg = 0
      c2g = 0
    else : print "model not implemented"
    print (s,float(float(t)/10.0))
    if 1>0  : #neventsTH >nev/0.25  and neventsTH <nev/0.01 : # neventsTH >nev/0.15 and neventsTH <nev/0.01 : 
      print "pottential to probe"
      os.system('python recastHH.py --kl '+str(float(s))+' --kt '+str(kt)+' --c2 '+str(c2)+' --cg '+str(cg)+' --c2g '+str(c2g)+' --exp '+exp+' --out '+name+' --bunch '+str(options.bunch)+\
                ' --par1 '+str(par1)+' --par2 '+str(par2))
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
