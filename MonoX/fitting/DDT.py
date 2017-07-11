#
import os
import ROOT
from ROOT import *
from array import array
import math
from math import *
import sys
import pdb

def ComputeDDT(name, point, nPtBins, nRhoBins, H):
	DDT = TH2F(name, "", nRhoBins, -6, -1.5, nPtBins, 200, 800)
	DDT.SetStats(0)
	nXb = H.GetXaxis().GetNbins()
	nYb = H.GetYaxis().GetNbins()
	for x in range(nXb):
		for y in range(nYb):
			proj = H.ProjectionZ("H3"+str(x)+str(y),x+1,x+1,y+1,y+1)
			print str(x+1) + "," + str(y+1) + ":    "+ str(proj.Integral())
			p = array('d', [point])
			q = array('d', [0.0]*len(p))
			proj.GetQuantiles(len(p), q, p)
			DDT.SetBinContent( x+1, y+1, q[0] );
	return DDT

def DisplayDDT(DDT, Title, SaveName):
	C = TCanvas("TempCanvas", "Title", 575, 500)
	C.cd()
	DDT.SetStats(0)
	DDT.GetXaxis().SetTitle("jet #rho")
	DDT.GetYaxis().SetTitle("jet p_{T}")
	DDT.Draw("COLZ")
	C.Print("MAP_"+SaveName+".gif")



#Bkgs =["/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/TTbar.root", "/data/t3home000/jorgem/lpc/mcremone/panda/v_8026_0_4/flat/control/WJets.root","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/ZJets.root"]
#Bkgs_tags=[("TTbar","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/TTbar.root"),("WJets","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/WJets.root"),("ZJets","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/ZJets.root")]
#Bkgs_tags=[("TTbar","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/TTbar.root"),("WJets","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/WJets.root"),("ZJets","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control/ZJets.root")]
Bkgs_tags=[("Diboson","/data/t3home000/mcremone/lpc/jorgem/panda/v_8026_0_4/flat/control///fitting/fittingForest_test.root")]
H3={}


for bks,B in Bkgs_tags:

	H3[bks] = TH3F("H3_%s"%(bks), "H3_%s"%(bks), 9, -6, -1.5, 12, 200, 800, 500, 0, 0.5)
	
for bks,B in Bkgs_tags:

	print 'Starting with '+bks+'-------------------------- :)'

	H3[bks].SetStats(0)
	F = TFile(B)
	if "test" in B:
		tree = "Diboson_test"
	T = F.Get("%s"%tree)
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
                if(j % (1 * n/100) == 0):
                        sys.stdout.write("\r[" + "="*int(20*j/n) + " " + str(round(100.*j/n,0)) + "% done")
                        sys.stdout.flush()
		T.GetEntry(j)
		#weight = T.puWeight*T.scale1fb*T.kfactor*T.kfactorNLO
		weight = 1
		PT = T.fj1Pt
		preRHO = T.fj1MSD_corr*T.fj1MSD_corr/T.fj1Pt/T.fj1Pt
		if preRHO > 0. and T.fj1ECFN_1_2_10 != 0.:
			RHO = math.log(preRHO)
			jtN2b1sd_8 = T.fj1ECFN_2_3_10/math.pow(T.fj1ECFN_1_2_10,2.00)
			if PT > 200 and RHO < -1.5 and RHO > -6.0 and jtN2b1sd_8 > 0.:
				H3[bks].Fill(RHO, PT, jtN2b1sd_8, weight)
DDT_5by6={}
DDT_5by3={}
Fout = TFile("DDTs.root", "recreate")
Fout.cd()
for key in H3:
	DDT_5by6[key]=ComputeDDT('DDT_5by6_%s'%(key), 0.05, 12, 9, H3[key])
	DDT_5by6[key].Write()
	DDT_5by3[key]=ComputeDDT('DDT_5by3_%s'%(key), 0.05, 3, 9, H3[key])
	DDT_5by3[key].Write()
Fout.Close()

#DDT_5by6 = ComputeDDT("DDT_5by6", 0.05, 12, 9, H3)
#DisplayDDT(DDT_5by6, "DDT vals at 5%", "DDT_5by6")
#DDT_5by3 = ComputeDDT("DDT_5by3", 0.05, 3, 9, H3)
#DisplayDDT(DDT_5by3, "DDT vals at 5%", "DDT_5by3")

