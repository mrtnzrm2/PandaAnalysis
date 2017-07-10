from PandaCore.Tools.Misc import *
from re import sub

triggers = {
    'met':'(trigger&1)!=0',
    'ele':'(trigger&2)!=0',
    'pho':'(trigger&4)!=0',
}

metFilter='metFilter==1 && egmFilter==1'
topTagSF = '1'
presel = 'nFatjet==1 && fj1Pt>200 && nTau==0 && Sum$(jetPt>30 && jetIso)<2 && fj1MSD_corr>100 && fj1MSD_corr<150'

cuts = {
    'signal'             : tAND(metFilter,tAND(presel,'nLooseLep==0 && nLooseElectron==0 && nLoosePhoton==0 && pfmet>250 && dphipfmet>0.4 && isojetNBtags==0')), 
    'singlemuon'         : tAND(metFilter,tAND(presel,'nLoosePhoton==0 && nTau==0 && nLooseLep==1 && looseLep1IsTight==1 && abs(looseLep1PdgId)==13 && pfUWmag>250 && dphipfUW>0.4 && mT<160')),
    'singleelectron'     : tAND(metFilter,tAND(presel,'nLoosePhoton==0 && nTau==0 && nLooseLep==1 && looseLep1IsTight==1 && looseLep1IsHLTSafe==1 && abs(looseLep1PdgId)==11 && pfmet>50 && pfUWmag>250 && dphipfUW>0.4 && mT<160')),
    'dimuon'            : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.4 && nLooseElectron==0 && nLoosePhoton==0 && nTau==0 && nLooseMuon==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && isojetNBtags==0')),
    'dielectron'        : tAND(metFilter,tAND(presel,'pfUZmag>250 && dphipfUZ>0.4 && nLooseMuon==0 && nLoosePhoton==0 && nTau==0 && nLooseElectron==2 && nTightLep>0 && 60<diLepMass && diLepMass<120 && isojetNBtags==0')),
#    'qcd'               : tAND(metFilter,tAND(removeCut(removeCut(presel,'top_ecf_bdt'),'fj1MSD'),'pfmet>250 && dphipfmet<0.1 && nLooseLep==0 && nLoosePhoton==0 && nTau==0 && fabs(calomet-pfmet)/pfmet<0.5')), 
}
for r in ['singlemuon','singleelectron']:
	cuts[r+'w'] = tAND(cuts[r],'isojetNBtags==0')
	cuts[r+'top'] = tAND(cuts[r],'isojetNBtags==1')


weights = {
  'signal'         : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_metTrig*sf_btag0',
  'top'            : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag1',
  'w'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag0',
  'notag'          : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV',
  'z'              : '%f*sf_pu*sf_tt*normalizedWeight*sf_lepID*sf_lepIso*sf_lepTrack*sf_ewkV*sf_qcdV*sf_btag0',
#  'photon'         : '%f*sf_pu*normalizedWeight*sf_ewkV*sf_qcdV*sf_pho*sf_phoTrig *sf_qcdV2j*sf_btag0', # add the additional 2-jet kfactor
}
weights['qcd'] = weights['signal']

for x in ['singlemuontop','singleelectrontop']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['top'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['top'],'sf_metTrig')
for x in ['singlemuonw','singleelectronw']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['w'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['w'],'sf_metTrig')
for x in ['dimuon','dielectron']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['z'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['z'],'sf_metTrig')
for x in ['singlemuon','singleelectron']:
	if 'electron' in x:
	  weights[x] = tTIMES(weights['notag'],'sf_eleTrig')
	else:
	  weights[x] = tTIMES(weights['notag'],'sf_metTrig')

for r in ['signal','top','w','singlemuontop','singleelectrontop','singlemuonw','singleelectronw','dimuon','dielectron','z']:
  for shift in ['BUp','BDown','MUp','MDown']:
    for cent in ['sf_btag']:
      weights[r+'_'+cent+shift] = sub(cent+'0',cent+'0'+shift,sub(cent+'1',cent+'1'+shift,weights[r]))
