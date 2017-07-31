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

basedir = getenv('SKIM_MONOHIGGS_BOOSTED_FLATDIR')+'/%s'%subfolder
lumi = 35900

def f(x):
    return basedir + x + '.root'

def transfile(process):
	name = subfolder.split("/")[0]
	trans = r.TFile(basedir+'/DDT_%s.root'%name)
	corrF56 = trans.Get("DDT_5by6_%s"%process)
	corrF53 = trans.Get("DDT_5by3_%s"%process)
	trans.Close()

	h = [corrF56,corrF53]
	return h
def addN2DDT(process):
	name = subfolder.split("/")[0]
	trans = r.TFile(basedir+'/DDT_%s.root'%name)
	f56 = trans.Get("DDT_5by6_%s"%process)
	f53 = trans.Get("DDT_5by3_%s"%process)
	fi = r.TFile(f(process),"update")
	fi.cd()
	t = fi.Get("events")
	#print t
	ndd56 = array('f', [0.0])
	t.Branch("n2ddt56", ndd56, 'n2ddt56'+'/F')
	ndd53 = array('f', [0.0])
	t.Branch("n2ddt53", ndd53, 'n2ddt53'+'/F')
	n = t.GetEntries()
	for j in range(0,n):
		t.GetEntry(j)
                pt = t.fj1Pt
		mass = t.fj1MSD_corr
		if t.fj1ECFN_1_2_10 == 0: continue
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
elif 'zee'  or 'zmn' in region:
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
		addN2DDT(process)

factory = forest.RegionFactory(name = region if not(is_test) else 'test',
                               cut = sel.cuts[region],
                               variables = vmap, 
                               mc_variables = mc_vmap, 
                               mc_weights = weights)

#Process and creation of new ntuples process
if is_test:
    factory.add_process(f('Diboson'),'Diboson')

#photon CR
elif region=='pho':
    factory.add_process(f('GJets'),'Pho')
    factory.add_process(f('SinglePhoton'),'Data',is_data=True)
    #factory.add_process(f('SinglePhoton'),'QCD',is_data=True,
    #                    extra_weights='sf_phoPurity',extra_cut=sel.triggers['pho'])

#photon fail CR
elif region=='pho_fail':
    factory.add_process(f('GJets'),'Pho')
    factory.add_process(f('SinglePhoton'),'Data',is_data=True,extra_cut=sel.phoTrigger)

    #factory.add_process(f('SinglePhoton'),'QCD',is_data=True,
    #                    extra_weights='sf_phoPurity',extra_cut=sel.triggers['pho'])

elif out_region not in ['signal_scalar','signal_vector','signal_thq','signal_stdm','signal']:

    factory.add_process(f('ZtoNuNu'),'Zvv')
    factory.add_process(f('ZJets'),'Zll')
    factory.add_process(f('WJets'),'Wlv')
    factory.add_process(f('SingleTop'),'ST')
    factory.add_process(f('Diboson'),'Diboson')
    factory.add_process(f('QCD'),'QCD')

    if ('zee' or 'te' or 'wen') == region:
        factory.add_process(f('SingleElectron'),'Data',is_data=True,extra_cut=sel.eleTrigger)
    	factory.add_process(f('TTbar'),'ttbar')


    if ('zee_fail' or 'te_fail' or 'wen_fail') == region:
        factory.add_process(f('SingleElectron'),'Data',is_data=True,extra_cut=sel.eleTrigger)


    if ('zmm' or 'tm' or 'tm_fail' or 'wmn') ==  region:
        factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.metTrigger)
    	factory.add_process(f('TTbar'),'ttbar')


    if ('zmm_fail' or 'wmn_fail')== region:
        factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.metTrigger)
	
elif out_region in ['signal_scalar','signal_vector','signal_thq','signal_stdm','signal']:
    factory.add_process(f('MET'),'Data',is_data=True,extra_cut=sel.metTrigger)
    factory.add_process(f('ZtoNuNu'),'Zvv')
    factory.add_process(f('TTbar'),'ttbar'[0])
    factory.add_process(f('ZJets'),'Zll')
    factory.add_process(f('WJets'),'Wlv')
    factory.add_process(f('SingleTop'),'ST')
    factory.add_process(f('Diboson'),'Diboson')
    factory.add_process(f('QCD'),'QCD')
    factory.add_process(f('ZpA0-600-300'),'ZpA0_600')
    factory.add_process(f('ZpA0-800-300'),'ZpA0_800')
    factory.add_process(f('ZpA0-1000-300'),'ZpA0_1000')
    factory.add_process(f('ZpA0-1200-300'),'ZpA0_1200')
    factory.add_process(f('ZpA0-1400-300'),'ZpA0_1400')
    factory.add_process(f('ZpA0-1700-300'),'ZpA0_1700')
    factory.add_process(f('ZpA0-2000-300'),'ZpA0_2000')
    factory.add_process(f('ZpA0-2500-300'),'ZpA0_2500')
    factory.add_process(f('ZpBaryonic-10-1'),'BarZp_10_1')
    factory.add_process(f('ZpBaryonic-10-10'),'BarZp_10_10')
    factory.add_process(f('ZpBaryonic-10-50'),'BarZp_10_50')
    factory.add_process(f('ZpBaryonic-10-150'),'BarZp_10_150')
    factory.add_process(f('ZpBaryonic-10-500'),'BarZp_10_500')
    factory.add_process(f('ZpBaryonic-15-10'),'BarZp_15_10')
    factory.add_process(f('ZpBaryonic-20-1'),'BarZp_20_1')
    factory.add_process(f('ZpBaryonic-50-1'),'BarZp_50_1')
    factory.add_process(f('ZpBaryonic-50-10'),'BarZp_50_10')
    factory.add_process(f('ZpBaryonic-50-50'),'BarZp_50_50')
    factory.add_process(f('ZpBaryonic-95-50'),'BarZp_95_50')
    factory.add_process(f('ZpBaryonic-100-1'),'BarZp_100_1')
    factory.add_process(f('ZpBaryonic-100-10'),'BarZp_100_10')
    factory.add_process(f('ZpBaryonic-200-1'),'BarZp_200_1')
    factory.add_process(f('ZpBaryonic-200-50'),'BarZp_200_50')
    factory.add_process(f('ZpBaryonic-200-150'),'BarZp_200_150')
    factory.add_process(f('ZpBaryonic-295-150'),'BarZp_295_150')
    factory.add_process(f('ZpBaryonic-300-1'),'BarZp_300_1')
    factory.add_process(f('ZpBaryonic-300-50'),'BarZp_300_50')
    factory.add_process(f('ZpBaryonic-500-1'),'BarZp_500_1')
    factory.add_process(f('ZpBaryonic-500-150'),'BarZp_500_150')
    factory.add_process(f('ZpBaryonic-995-500'),'BarZp_995_500')
    factory.add_process(f('ZpBaryonic-1000-1'),'BarZp_1000_1')
    factory.add_process(f('ZpBaryonic-1000-150'),'BarZp_1000_150')
    factory.add_process(f('ZpBaryonic-1000-1000'),'BarZp_1000_1000')
    factory.add_process(f('ZpBaryonic-1995-1000'),'BarZp_1995_1000')
    factory.add_process(f('ZpBaryonic-2000-1'),'BarZp_2000_1')
    factory.add_process(f('ZpBaryonic-10000-1'),'BarZp_10000_1')
    factory.add_process(f('ZpBaryonic-10000-10'),'BarZp_10000_10')
    factory.add_process(f('ZpBaryonic-10000-50'),'BarZp_10000_50')
    factory.add_process(f('ZpBaryonic-10000-150'),'BarZp_10000_150')
    factory.add_process(f('ZpBaryonic-10000-500'),'BarZp_10000_500')
'''
elif out_region=='signal_vector':
    signal_files = glob(basedir+'/Vector*root')
    if couplings:
        out_region += '_'+couplings
    for f in signal_files:
        fname = f.split('/')[-1].replace('.root','')
        signame = fname
        replacements = {
            'Vector_MonoTop_NLO_Mphi-':'',
            '_gSM-0p25_gDM-1p0_13TeV-madgraph':'',
            '_Mchi-':'_',
        }
        for k,v in replacements.iteritems():
            signame = signame.replace(k,v)
        factory.add_process(f,signame)
elif out_region=='signal_scalar':
    signal_files = glob(basedir+'/Scalar*root')
    if couplings:
        out_region += '_'+couplings
    for f in signal_files:
        fname = f.split('/')[-1].replace('.root','')
        signame = fname
        replacements = {
            'Scalar_MonoTop_LO_Mphi-':'',
            '_13TeV-madgraph':'',
            '_Mchi-':'_',
        }
        for k,v in replacements.iteritems():
            signame = signame.replace(k,v)
        factory.add_process(f,'scalar_'+signame)
elif out_region=='signal_thq':
    factory.add_process(f('thq'),'thq')
elif out_region=='signal_stdm':
    for m in [300,500,1000]:
        factory.add_process(f('ST_tch_DM-scalar_LO-%i_1-13_TeV'%m),'stdm_%i'%m)
'''
forestDir = '/data/t3home000/mcremone/lpc/jorgem/skim/monohiggs_boosted/'
os.system('mkdir -p %s/%s'%(forestDir,'fittingForest'))
factory.run(forestDir+'/fittingForest/fittingForest_%s.root'%out_region)

