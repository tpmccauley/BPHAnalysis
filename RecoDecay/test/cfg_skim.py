import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring(
#
### use this to access the nearest copy of the input file, querying the catalog
#
    '/store/group/lpcmuon/dimuon/MuOnia/Run2015D-Onia2MuMuPAT-v4/151209_052543/0000/Onia2MuMuPAT_1.root'
### use this to access the input file if by any reason you want to specify 
### the data server
#    'root://xrootd-cms.infn.it//store/group/lpcmuon/dimuon/MuOnia/Run2015D-Onia2MuMuPAT-v4/151209_052543/0000/Onia2MuMuPAT_1.root'
#
### use this to access an input file locally available
#    'file:/...complete_file_path.../Onia2MuMuPAT_1.root'
))

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.testBPHRecoDecay = cms.EDAnalyzer('TestBPHRecoDecay',
    gpCandsLabel = cms.string('cleanPatTrackCands::PAT'),
    ccCandsLabel = cms.string('onia2MuMuPAT::PAT'),
    outDump = cms.string('dump_skim.txt'),
    outHist = cms.string('hist_skim.root')
)

process.p = cms.Path(
    process.testBPHRecoDecay
)

