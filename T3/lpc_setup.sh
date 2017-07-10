#!/bin/bash                                                                                                                                                                                                 
export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/
export SUBMIT_NAME="v_8026_0_4"

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
#export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170316.cfg"
#export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20170501_ttdm.cfg"
export PANDA_CFG="http://t3serv001.mit.edu/~mcremone/histcatalog/test.cfg"
export PANDA_FLATDIR="/uscms_data/d3/matteoc/panda/"${SUBMIT_NAME}"/flat/"
#export PANDA_FLATDIR="/data/t3home000/mcremone/lpc/jorgem/panda/"${SUBMIT_NAME}"/flat/control/"                                                                                                            
mkdir -p $PANDA_FLATDIR

#export SUBMIT_TMPL="skim_monojet_tmpl.py" ####
export SUBMIT_TMPL="skim_vbf_tmpl.py"
                                                                                                                                            
export SUBMIT_WORKDIR="/uscms_data/d3/matteoc/condor/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/uscms_data/d3/matteoc/condor/"${SUBMIT_NAME}"/logs/"
export SUBMIT_OUTDIR="/store/user/matteoc/panda/"${SUBMIT_NAME}"/batch/"

export SKIM_CFGDIR="/uscms_data/d3/matteoc/skim/configs"

export SKIM_MONOJET_FLATDIR="/uscms_data/d3/matteoc/skim/"${SUBMIT_NAME}"/monojet/"
export SKIM_MONOHIGGS_FLATDIR="/uscms_data/d3/matteoc/skim/"${SUBMIT_NAME}"/monohiggs_boosted/"
export SKIM_MONOHIGGS_RESOLVED_FLATDIR="/uscms_data/d3/matteoc/skim/"${SUBMIT_NAME}"/monohiggs_resolved/"

mkdir -p $SUBMIT_WORKDIR $SUBMIT_LOGDIR $SKIM_MONOJET_FLATDIR $SKIM_MONOHIGGS_FLATDIR $SKIM_MONOHIGGS_RESOLVED_FLATDIR
eosmkdir -p $SUBMIT_OUTDIR/locks/

rm $SKIM_MONOJET_FLATDIR/*.sh
rm $SKIM_MONOHIGGS_FLATDIR/*.sh
rm $SKIM_MONOHIGGS_RESOLVED_FLATDIR/*.sh

ln -s $SKIM_CFGDIR/runSkim.sh  $SKIM_MONOJET_FLATDIR
ln -s $SKIM_CFGDIR/runSkimAll.sh  $SKIM_MONOJET_FLATDIR

ln -s $SKIM_CFGDIR/runSkim.sh  $SKIM_MONOHIGGS_FLATDIR
ln -s $SKIM_CFGDIR/runSkimAll.sh $SKIM_MONOHIGGS_FLATDIR

ln -s $SKIM_CFGDIR/runSkim.sh  $SKIM_MONOHIGGS_RESOLVED_FLATDIR
ln -s $SKIM_CFGDIR/runSkimAll.sh $SKIM_MONOHIGGS_RESOLVED_FLATDIR

#private production                                                                                                                                                                                         
export PRIVATE_LOGDIR="${HOME}/cms/logs/monotop_private_panda/"
export PRIVATE_PRODDIR="${HOME}/cms/hist/monotop_private_pandatree/"
export PRIVATE_CFGDIR="${HOME}/cms/condor/monotop_private_panda/"

# fitting                                                                                                                                                                                                   
#export PANDA_FIT=/data/t3serv014/snarayan/CMSSW_7_4_7/
#export PANDA_XSECS=/home/snarayan/cms/cmssw/analysis/MonoTop_Xsec/
#export PANDA_FITTING=${PANDA_FLATDIR}/fitting/
#mkdir -p $PANDA_FITTING/scans/ $PANDA_FITTING/logs/

