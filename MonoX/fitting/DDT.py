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
	DDT = TH2F(name, "", nRhoBins, -6, -1, nPtBins, 175, 825)
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

H3 = TH3F("H3", "", 10, -6, -1, 13, 175, 825, 500, 0, 0.5)
H3.SetStats(0)
Bkgs =["QCDmvaEVwithExtv3.root", "GJetsmvaEVv3.root"]
for B in Bkgs:
	F = TFile(B)
	T = F.Get("Events")
	n = T.GetEntries()
	for j in range(0, n): # Here is where we loop over all events.
		T.GetEntry(j)
		weight = T.puWeight*T.scale1fb*T.kfactor*T.kfactorNLO
		if T.vpho0_pt>175 and T.vpho0_mva > 0.5 and T.AK8Puppijet0_pt>175: # photon requirement
			PT = T.AK8Puppijet0_pt
			preRHO = T.AK8Puppijet0_msd*T.AK8Puppijet0_msd/T.AK8Puppijet0_pt/T.AK8Puppijet0_pt
			if preRHO > 0.:
				RHO = math.log(preRHO)
				if RHO < -1 and RHO > -6 and T.AK8Puppijet0_N2sdb1 > 0.:
					t3pho = TVector3()
					t3pho.SetPtEtaPhi(T.vpho0_pt, T.vpho0_eta, T.vpho0_phi)
					t3jet = TVector3()
					t3jet.SetPtEtaPhi(T.AK8Puppijet0_pt, T.AK8Puppijet0_eta, T.AK8Puppijet0_phi)
					if t3pho.DeltaR(t3jet) > 2.2:
						H3.Fill(RHO, PT, T.AK8Puppijet0_N2sdb1, weight)

DDT_photon = ComputeDDT("DDT_photon", 0.05, 13, 10, H3)
DisplayDDT(DDT_photon, "DDT vals at 5%", "DDT_photon")

Fout = TFile("PhotonDDTs.root", "recreate")
Fout.cd()
DDT_photon.Write()
Fout.Close()
