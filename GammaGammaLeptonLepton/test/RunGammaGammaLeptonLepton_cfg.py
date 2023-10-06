import FWCore.ParameterSet.Config as cms

process = cms.Process("ggll")

runOnMC = False
DORERECO = True # Set this to run on-the-fly re-reco of protons with new conditions

#########################
#    General options    #
#########################

process.load("FWCore.MessageService.MessageLogger_cfi")
process.options   = cms.untracked.PSet(
    #wantSummary = cms.untracked.bool(True),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    allowUnscheduled = cms.untracked.bool(True),
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

#########################
#      Input files      #
#########################

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2018D/DoubleMuon/AOD/15Feb2022_UL2018-v1/50000/B2D2E689-C44D-4B43-ABF9-2B5CFA3B2F5D.root',
        '/store/data/Run2018D/DoubleMuon/AOD/15Feb2022_UL2018-v1/50000/B1C2608E-C8C4-084D-BF8D-3D89D9C9B929.root',
        '/store/data/Run2018D/DoubleMuon/AOD/15Feb2022_UL2018-v1/50000/B2D2E689-C44D-4B43-ABF9-2B5CFA3B2F5D.root',
        '/store/data/Run2018D/DoubleMuon/AOD/15Feb2022_UL2018-v1/50000/BCF36728-F59E-2544-B0DC-EFEC35375552.root',
        '/store/data/Run2018D/DoubleMuon/AOD/15Feb2022_UL2018-v1/50000/BECEABB6-78A7-B741-B166-0DECEB5CD7FD.root',
        '/store/data/Run2018D/DoubleMuon/AOD/15Feb2022_UL2018-v1/50000/BFA0B870-B6E3-E44F-B0DD-124B5322C521.root'

    ),
    #firstEvent = cms.untracked.uint32(0)
)

#########################
#        Triggers       #
#########################

process.load("TaggedProtonDileptonValidation.GammaGammaLeptonLepton.HLTFilter_cfi")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilter.HLTPaths = cms.vstring(
    'HLT_DoubleMu43NoFiltersNoVtx_*',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_*',
    'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_*'
)

#########################
#      Preskimming      #
#########################
process.load("Configuration.StandardSequences.GeometryDB_cff") 
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

#########################
#     Proton RECO       #
#########################
process.load("RecoPPS.Configuration.recoCTPPS_cff")

if DORERECO == True:
    process.ctppsLocalTrackLiteProducer.tagPixelTrack = cms.InputTag("ctppsPixelLocalTracks","","ggll")
    process.ctppsLocalTrackLiteProducer.tagDiamondTrack = cms.InputTag("ctppsDiamondLocalTracks","","ggll")
    process.ctppsProtons.tagLocalTrackLite = cms.InputTag("ctppsLocalTrackLiteProducer","","ggll")
    process.ctppsLocalTrackLiteProducer.includePixels = cms.bool(True)
    process.ctppsLocalTrackLiteProducer.includeDiamonds = cms.bool(True)
    process.ctppsProtons.doSingleRPReconstruction = cms.bool(True)
    process.ctppsProtons.doMultiRPReconstruction = cms.bool(True)

####################################################################
# Change this part to pick up new Alignment/Optics conditions 
# for re-reconsturcing protons
####################################################################

#process.load("RecoCTPPS.ProtonReconstruction.year_2017_OF.ctppsProtonReconstructionOF_cfi")
#process.load("RecoCTPPS.ProtonReconstruction.year_2018_OFDB.ctppsProtonReconstructionOFDB_cfi")
# conditions DB for 2018                                                                                                                    
                         
#from CondCore.CondDB.CondDB_cfi import *
#
#CondDB.connect = 'frontier://FrontierProd/CMS_CONDITIONS'
#
#process.PoolDBESSource2 = cms.ESSource("PoolDBESSource",
#                                       CondDB,
#                                       DumpStat = cms.untracked.bool(False),
#                                       toGet = cms.VPSet(cms.PSet(
#            record = cms.string('LHCInfoRcd'),
#            #tag = cms.string("LHCInfoTest_prompt_v3")  
#            tag = cms.string("LHCInfoEndFill_prompt_v1")
#            )),
#                                       )
#
## get optics from Wagners's DB tag
#from CondCore.CondDB.CondDB_cfi import *
#process.CondDBOptics = CondDB.clone( connect = 'frontier://FrontierProd/CMS_CONDITIONS' )
#process.PoolDBESSourceOptics = cms.ESSource("PoolDBESSource",
#    process.CondDBOptics,
#    DumpStat = cms.untracked.bool(False),
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('CTPPSOpticsRcd'),
#        tag = cms.string("PPSOpticalFunctions_offline_v1")
#    )),
#)
#
## get alignment from Clemencia's file
#from CondCore.CondDB.CondDB_cfi import *
#process.CondDBAlignment = CondDB.clone( connect = 'sqlite_file:/afs/cern.ch/user/c/cmora/public/CTPPSDB/AlignmentSQlite/CTPPSR#PRealAlignment_table_v26Apr.db' )
#process.PoolDBESSourceAlignment = cms.ESSource("PoolDBESSource",
#    process.CondDBAlignment,
#    #timetype = cms.untracked.string('runnumber'),
#    toGet = cms.VPSet(cms.PSet(
#        record = cms.string('RPRealAlignmentRecord'),
#        tag = cms.string('CTPPSRPAlignment_real_table_v26A19')
#    ))
#)
#
##JH - ESPrefer to get optical functions from CTPPSOpticalFunctionsESSource instead of global tag for now                     
# process.es_prefer_ppsOptics = cms.ESPrefer("PoolDBESSource","PoolDBESSourceOptics")


#########################
#       Analysis        #
#########################

process.load("TaggedProtonDileptonValidation.GammaGammaLeptonLepton.GammaGammaLL_cfi")

process.ggll_aod.triggersList = process.hltFilter.HLTPaths
process.ggll_aod.leptonsType = cms.string('Muon')
process.ggll_aod.runOnMC = cms.bool(runOnMC)
process.ggll_aod.fetchProtons = cms.bool(True)
process.ggll_aod.saveExtraTracks = cms.bool(False)
process.ggll_aod.year = cms.string('2018')
if DORERECO == True:
    # Use protons re-recoed with new conditions
    process.ggll_aod.ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP","ggll")
    process.ggll_aod.ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP","ggll")
else:
    # Use protons directly from the AOD
    process.ggll_aod.ppsRecoProtonSingleRPTag = cms.InputTag("ctppsProtons", "singleRP")
    process.ggll_aod.ppsRecoProtonMultiRPTag = cms.InputTag("ctppsProtons", "multiRP")
process.ggll_aod.lhcInfoLabel = cms.string('')
process.ggll_aod.stageL1Trigger = cms.uint32(2)

# prepare the output file
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('output.root'),
    closeFileFast = cms.untracked.bool(True)
)

if DORERECO == True:
    process.p = cms.Path(
        process.hltFilter*
        process.ctppsDiamondLocalTracks*
        process.ctppsPixelLocalReconstruction *
        process.ctppsLocalTrackLiteProducer *
        process.ctppsProtons * 
        process.ggll_aod
    )
else:
    process.p = cms.Path(
        process.hltFilter*
        process.ggll_aod
    )



#print process.dumpPython()
