#!/usr/bin/env python
from re import sub
from math import *
from array import array
from sys import argv,exit
from os import path,getenv
import os
import ROOT as r
from glob import glob
import argparse
parser = argparse.ArgumentParser(description='make forest')
parser.add_argument('--region',metavar='region',type=str,default=None)
parser.add_argument('--couplings',metavar='couplings',type=str,default=None)
parser.add_argument('--var',metavar='var',type=str,default=None)
parser.add_argument('--input',metavar='input',type=str,default='/data/t3home000/mcremone/lpc/jorgem/skim/v_8026_0_4/monohiggs_boosted/')
parser.add_argument('--output',metavar='output',type=str,default='/data/t3home000/mcremone/lpc/jorgem/skim/monohiggs_boosted/')
parser.add_argument('--cr',metavar='cr',type=str,default=None)
args = parser.parse_args()
couplings = args.couplings
nddt = args.var
if couplings=='nominal':
    couplings = None
out_region = args.region
region = out_region.split('_')[0]
if region=='test':
    is_test = True 
    region = 'signal'
else:
    is_test = False

argv=[]
import PandaAnalysis.Flat.fitting_forest as forest 
from PandaCore.Tools.Misc import *
import PandaCore.Tools.Functions # kinematics
import PandaAnalysis.MonoH.Selection_doubleb as sel
#import PandaAnalysis.Monotop.CombinedSelection as sel



toProcess = out_region

#Addtional selection

weights_ttbar = {'ttbar_mistag_SF':"1.03*"}
weights_signal = {'signal_eff_SF':"((fj1Pt<350)*0.91+(fj1Pt>350)*1)*"}
otherWeights = [weights_ttbar , weights_signal]

#selection of process
subfolder = ''
if toProcess=='signal':
  subfolder = 'sr/'
elif toProcess=='wmn':
  subfolder = 'cr_w_mu/'
elif toProcess=='wmn_fail':
  subfolder = 'cr_w_mu/'
elif toProcess=='wen':
  subfolder = 'cr_w_el/'
elif toProcess=='wen_fail':
  subfolder = 'cr_w_el/'
elif toProcess=='tm':
  subfolder = 'cr_ttbar_mu/'
elif toProcess=='tm_fail':
  subfolder = 'cr_ttbar_mu/'
elif toProcess=='te':
  subfolder = 'cr_ttbar_el/'
elif toProcess=='te_fail':
  subfolder = 'cr_ttbar_el/'
elif toProcess=='zmm':
  subfolder = 'cr_dimuon/'
elif toProcess=='zmm_fail':
  subfolder = 'cr_dimuon/'
elif toProcess=='zee':
  subfolder = 'cr_dielectron/'
elif toProcess=='zee_fail':
  subfolder = 'cr_dielectron/'
elif toProcess=='pho':
  subfolder = 'cr_gamma/'
elif toProcess=='pho_fail':
  subfolder = 'cr_gamma/'

#basedir = getenv('SKIM_MONOHIGGS_BOOSTED_FLATDIR')+'/%s'%subfolder
basedir = args.input
if args.cr:
	basedir = '%s/%s'%(args.input,subfolder)
lumi = 35900

def f(x):
    return basedir+'/' + x + '.root'

def addN2DDT(process):
	name = subfolder.split("/")[0]
	trans = r.TFile(basedir+'/DDT.root')	#basedir needs to point to a folder where the DDT mapping is present
	f56 = trans.Get("DDT_5by6_%s"%process)
	f53 = trans.Get("DDT_5by3_%s"%process)
	fi = r.TFile(f(process),"update")
	fi.cd()
	t = fi.Get("events")
	ndd56 = array('f', [0.0])
	t.Branch("n2ddt56", ndd56, 'n2ddt56'+'/F')
	ndd53 = array('f', [0.0])
	t.Branch("n2ddt53", ndd53, 'n2ddt53'+'/F')
	n = t.GetEntries()
	for j in range(0,n):
		t.GetEntry(j)
                pt = t.fj1Pt
		mass = t.fj1MSD_corr
		if t.fj1ECFN_1_2_10 == 0: 
			
			ndd56[0] = 999
			ndd53[0] = 999
			t.Fill()
			continue
		jtN2b1sd_8 = t.fj1ECFN_2_3_10/pow(t.fj1ECFN_1_2_10,2.00)
                if (jtN2b1sd_8 > 0. and mass*mass/pt/pt) > 0.0 and pt > 200 and pt < 800:

	                rho = log(mass*mass/pt/pt)
			rind6 = f56.GetXaxis().FindBin(rho)
                        pind6 = f56.GetYaxis().FindBin(pt)

                        rind3 = f53.GetXaxis().FindBin(rho)
                        pind3 = f53.GetYaxis().FindBin(pt)		
		
			if rho >  f56.GetXaxis().GetBinUpEdge( f56.GetXaxis().GetNbins() ) :
                        	rind6 = f56.GetXaxis().GetNbins()
                        if rho <  f56.GetXaxis().GetBinLowEdge( 1 ) :
                                rind6 = 1
                        if pt >  f56.GetYaxis().GetBinUpEdge( f56.GetYaxis().GetNbins() ) :
                                pind6 = f56.GetYaxis().GetNbins()
                        if pt < f56.GetYaxis().GetBinLowEdge( 1 ) :
                                pind6 = 1

                        if rho > f53.GetXaxis().GetBinUpEdge( f53.GetXaxis().GetNbins() ) :
                                rind3 = f53.GetXaxis().GetNbins()
                        if rho <  f53.GetXaxis().GetBinLowEdge( 1 ) :
                                rind3 = 1
                        if pt > f53.GetYaxis().GetBinUpEdge( f53.GetYaxis().GetNbins() ) :
                                pind3 = f53.GetYaxis().GetNbins()
                        if pt < f53.GetYaxis().GetBinLowEdge( 1 ) :
                                pind3 = 1

                        ndd56[0] = jtN2b1sd_8 - f56.GetBinContent(rind6,pind6)
                        ndd53[0] = jtN2b1sd_8 - f53.GetBinContent(rind3,pind3)
                        t.Fill()
		else:
			ndd56[0] = 999
                        ndd53[0] = 999
			t.Fill()

	fi.Write('',r.TObject.kWriteDelete)
	fi.Close()
	trans.Close()



def shift_btags(additional=None):
    shifted_weights = {}
    #if not any([x in region for x in ['signal','top','w']]):
    #    return shifted_weights 
    for shift in ['BUp','BDown','MUp','MDown']:
        for cent in ['sf_btag']:
            shiftedlabel = ''
            if 'sj' in cent:
                shiftedlabel += 'sj'
            if 'B' in shift:
                shiftedlabel += 'btag'
            else:
                shiftedlabel += 'mistag'
            if 'Up' in shift:
                shiftedlabel += 'Up'
            else:
                shiftedlabel += 'Down'
            weight = sel.weights[region+'_'+cent+shift]%lumi
            if additional:
                weight = tTIMES(weight,additional)
            shifted_weights[shiftedlabel] = weight
    return shifted_weights


#vmap definition
vmap = {}
mc_vmap = {'genBosonPt':'genBosonPt'}
if region in ['signal','test']:
    u,uphi, = ('puppimet','dphipuppimet')
elif 'pho' in region:
    u,uphi = ('puppiUAmag','dphipuppiUA')
elif 'wmn'or 'wen' or 'te' or 'tm' in region:
    u,uphi = ('puppiUWmag','dphipuppiUW')
elif 'zee'  or 'zmm' in region:
    u,uphi = ('puppiUZmag','dphipuppiUZ')
vmap['met'] = 'min(%s,999.9999)'%u 
vmap['ndd56'] = 'n2ddt56'
vmap['ndd53'] = 'n2ddt53'



#weights

weights = {'nominal' : sel.weights[region]%lumi}
if couplings:
    weights['nominal'] = tTIMES(weights['nominal'],couplings)
weights.update(shift_btags(couplings))

#Computing n2ddt variables in the ntuples stored inside the basedir directory
if nddt:
	for process in os.listdir(basedir):
		if 'DDT' in process: continue
		if ".root" in process: process = process.split(".")[0]
		print 'Starting process----',process
		addN2DDT(process)
		print ':) Added N2DDT --- Successfuly'


ps ={'ZtoNuNu':'Zvv','ZJets':'Zll','WJets':'Wlv','SingleTop':'ST','Diboson':'Diboson','QCD':'QCD','SingleElectron':'Data','TTbar':'ttbar','MET':'Data', 'GJets':'Pho' , 'SinglePhoton':'Data'}

forestDir = args.output
#if 'fail' in toProcess: os.system('mkdir -p %s/%s/%s'%(forestDir,'fittingForest',subfolder+'_fail'))
#else: os.system('mkdir -p %s/%s/%s'%(forestDir,'fittingForest',subfolder))
os.system('mkdir -p %s/%s/%s'%(forestDir,'fittingForest',subfolder))

for rfiles in os.listdir(basedir):
	if '.root' in rfiles: rfile = rfiles.split('.')[0]
	if rfile not in ps: continue
	if  rfile == 'SingleElectron' and ('zee' in region or 'te' in region or 'wen' in region):	
		print rfile,0
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[out_region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , ps[rfile] , is_data=True, extra_cut=sel.eleTrigger )
		if 'fail' in toProcess: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s_fail.root'%(subfolder,rfile))
		else: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s.root'%(subfolder,rfile))
	if rfile == 'MET' and ('zmm' in region or 'tm' in region or 'wmn' in region):
		print rfile,1
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[out_region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , ps[rfile] , is_data=True, extra_cut=sel.metTrigger )
		if 'fail' in toProcess: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s_fail.root'%(subfolder,rfile))
		else: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s.root'%(subfolder,rfile))
	if rfile == 'SinglePhoton' and  'pho' in region:
		print rfile,2
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[out_region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , ps[rfile] , is_data=True)
		if 'fail' in toProcess: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s_fail.root'%(subfolder,rfile))
                else: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s.root'%(subfolder,rfile))
	if ps[rfile] is not 'Data' and rfile is not ('SingleElectron' or 'SinglePhoton' or 'MET'):
		print rfile,3
		factory = forest.RegionFactory(name = region if not(is_test) else 'test',
        	                       cut = sel.cuts[out_region],
                	               variables = vmap, 
                        	       mc_variables = mc_vmap, 
                              	       mc_weights = weights)
		factory.add_process(f(rfile) , ps[rfile] )
		if 'fail' in toProcess: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s_fail.root'%(subfolder,rfile))
                else: factory.run(forestDir+'/fittingForest/%s/fittingForest_%s.root'%(subfolder,rfile))
	
pd = {'ZpA0-600-300':'ZpA0_600','ZpA0-800-300':'ZpA0_800','ZpA0-1000-300':'ZpA0_1000','ZpA0-1200-300':'ZpA0_1200','ZpA0-1400-300':'ZpA0_1400','ZpA0-1700-300':'ZpA0_1700','ZpA0-2000-300':'ZpA0_2000','ZpA0-2500-300':'ZpA0_2500','ZpBaryonic-10-1':'BarZp_10_1','ZpBaryonic-10-10':'BarZp_10_10','ZpBaryonic-10-50':'BarZp_10_50','ZpBaryonic-10-150':'BarZp_10_150','ZpBaryonic-10-500':'BarZp_10_500','ZpBaryonic-15-10':'BarZp_15_10','ZpBaryonic-20-1':'BarZp_20_1','ZpBaryonic-50-1':'BarZp_50_1','ZpBaryonic-50-10':'BarZp_50_10','ZpBaryonic-50-50':'BarZp_50_50','ZpBaryonic-95-50':'BarZp_95_50','ZpBaryonic-100-1':'BarZp_100_1','ZpBaryonic-100-10':'BarZp_100_10','ZpBaryonic-200-1':'BarZp_200_1','ZpBaryonic-200-50':'BarZp_200_50','ZpBaryonic-200-150':'BarZp_200_150','ZpBaryonic-295-150':'BarZp_295_150','ZpBaryonic-300-1':'BarZp_300_1','ZpBaryonic-300-50':'BarZp_300_50','ZpBaryonic-500-1':'BarZp_500_1','ZpBaryonic-500-150':'BarZp_500_150','ZpBaryonic-995-500':'BarZp_995_500','ZpBaryonic-1000-1':'BarZp_1000_1','ZpBaryonic-1000-150':'BarZp_1000_150','ZpBaryonic-1000-1000':'BarZp_1000_1000','ZpBaryonic-1995-1000':'BarZp_1995_1000','ZpBaryonic-2000-1':'BarZp_2000_1','ZpBaryonic-2000-500':'BarZp_2000_500','ZpBaryonic-10000-1':'BarZp_10000_1','ZpBaryonic-10000-150':'BarZp_10000_150','ZpBaryonic-10000-500':'BarZp_10000_500','ZpBaryonic-10000-10':'BarZp_10000_10','ZpBaryonic-10000-50':'BarZp_10000_50'}
	
if out_region in ['signal_scalar','signal_vector','signal_thq','signal_stdm','signal','sr']:

  for rfiles in os.listdir(basedir):
        if '.root' in rfiles: rfile = rfiles.split('.')[0]
        if rfile not in pd: continue
	if rfile is 'SingleElectron' or rfile is 'SinglePhoton': continue
        if rfile == 'MET' and ('zmm' in region or 'tm' in region or 'wmn' in region or 'sr' in region):
                #print rfile
                factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                                       cut = sel.cuts[out_region],
                                       variables = vmap,
                                       mc_variables = mc_vmap,
                                       mc_weights = weights)
                factory.add_process(f(rfile) , pd[rfile] , is_data=True, extra_cut=sel.metTrigger )
        	factory.run(forestDir+'/fittingForest/%s/fittingForest_%s.root'%(subfolder,rfile))
        else:
                #print rfile
                factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                                       cut = sel.cuts[out_region],
                                       variables = vmap,
                                       mc_variables = mc_vmap,
                                       mc_weights = weights)
                factory.add_process(f(rfile) , pd[rfile] )
 	        factory.run(forestDir+'/fittingForest/%s/fittingForest_%s.root'%(subfolder,rfile))

