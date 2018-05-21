import FWCore.ParameterSet.Config as cms

process = cms.Process("bphAnalysis")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#filename = open('/lustre/cmswork/ronchese/pat_ntu/cmssw910/chkMini02/src/PDAnalysis/EDM/prod/datafiles_Vmini/mc/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen_mini910/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen_mini910_020.list', 'r')
#fileList = cms.untracked.vstring( filename.readlines() )

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

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.options.allowUnscheduled = cms.untracked.bool(True)

process.source = cms.Source ("PoolSource", 
# Charmonium
    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2018A/Charmonium/MINIAOD/PromptReco-v2/000/316/239/00000/08BFAB4F-1359-E811-B01F-FA163E5285EC.root')
# DoubleMuon
#    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2018A/DoubleMuon/MINIAOD/PromptReco-v2/000/316/239/00000/3AC438F4-1859-E811-8996-FA163ECC1FE3.root')
# DoubleMuonLowMass
#    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2018A/DoubleMuonLowMass/MINIAOD/PromptReco-v2/000/316/239/00000/12F96450-1A59-E811-B49C-FA163E071950.root')
# MuOnia
#    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2018A/MuOnia/MINIAOD/PromptReco-v2/000/316/239/00000/027C6081-2A59-E811-A905-FA163E86C8EE.root')
#    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2017F/Charmonium/MINIAOD/17Nov2017-v1/00000/006425F0-6DED-E711-850C-0025904C66E8.root')
#    fileNames = cms.untracked.vstring('file:/home/ronchese/data/store/data/Run2017C/Charmonium/MINIAOD/PromptReco-v1/000/299/368/00000/64EE2B4F-866D-E711-A9DF-02163E01A6CE.root')
#    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2017B/Charmonium/MINIAOD/PromptReco-v2/000/298/853/00000/0C3845CC-8A68-E711-AAC9-02163E0118E2.root')
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9', '')

### vtxTagInfo
process.load('RecoBTag/SoftLepton/softLepton_cff')
#process.load('RecoBTag/SoftLepton/softPFMuonTagInfos_cfi')
#process.load('RecoBTag/SoftLepton/softPFElectronTagInfos_cfi')
process.load('RecoBTag/SecondaryVertex/pfInclusiveSecondaryVertexFinderTagInfos_cfi')
process.load('RecoBTag/ImpactParameter/pfImpactParameterTagInfos_cfi')

#process.load('RecoBTag/SecondaryVertex/secondaryVertexTagInfos_cfi')

process.softPFMuonsTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.softPFMuonsTagInfos.jets = cms.InputTag("slimmedJets")
process.softPFMuonsTagInfos.muons = cms.InputTag("slimmedMuons")

process.softPFElectronsTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.softPFElectronsTagInfos.jets = cms.InputTag("slimmedJets")
process.softPFElectronsTagInfos.electrons = cms.InputTag("slimmedElectrons")

process.pfImpactParameterTagInfos.primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfImpactParameterTagInfos.jets = cms.InputTag("slimmedJets")
process.pfImpactParameterTagInfos.candidates = cms.InputTag("packedPFCandidates")

process.tagInfoProd = cms.Sequence(
    process.softPFMuonsTagInfos
    + process.softPFElectronsTagInfos
    + process.pfImpactParameterTagInfos
    * process.pfSecondaryVertexTagInfos
    )

### vtxTagInfo end

#from BPHAnalysis.SpecificDecay.recoSelectForWrite_cfi import recoSelect
from BPHAnalysis.SpecificDecay.recoSelectForWrite_cfi import recoSelect

process.bphWriteSpecificDecay = cms.EDProducer('BPHWriteSpecificDecay',
    pVertexLabel = cms.string('offlineSlimmedPrimaryVertices'),
    pcCandsLabel = cms.string('packedPFCandidates::RECO'),
    patMuonLabel = cms.string('slimmedMuons'),
    kSCandsLabel = cms.string('slimmedKshortVertices::RECO'),
    lSCandsLabel = cms.string('slimmedLambdaVertices::RECO'),
    oniaName = cms.string('oniaFitted'),
    sdName   = cms.string('kx0Cand'),
    ssName   = cms.string('phiCand'),
    buName   = cms.string('buFitted'),
    bdName   = cms.string('bdFitted'),
    bsName   = cms.string('bsFitted'),
    k0Name   = cms.string('k0Fitted'),
    l0Name   = cms.string('l0Fitted'),
    b0Name   = cms.string('b0Fitted'),
    lbName   = cms.string('lbFitted'),
    writeVertex   = cms.bool( True ),
    writeMomentum = cms.bool( True ),
    recoSelect = cms.VPSet(recoSelect)
)

process.TFileService = cms.Service('TFileService',
  fileName = cms.string('his.root'),
  closeFileFast = cms.untracked.bool(True)
)

process.bphHistoSpecificDecay = cms.EDAnalyzer('BPHHistoSpecificDecay',
    trigResultsLabel  = cms.string('TriggerResults::HLT'), 
    oniaCandsLabel = cms.string('bphWriteSpecificDecay:oniaFitted:bphAnalysis'),
    sdCandsLabel = cms.string('bphWriteSpecificDecay:kx0Cand:bphAnalysis'),
    ssCandsLabel = cms.string('bphWriteSpecificDecay:phiCand:bphAnalysis'),
    buCandsLabel = cms.string('bphWriteSpecificDecay:buFitted:bphAnalysis'),
    bdCandsLabel = cms.string('bphWriteSpecificDecay:bdFitted:bphAnalysis'),
    bsCandsLabel = cms.string('bphWriteSpecificDecay:bsFitted:bphAnalysis'),
    k0CandsLabel = cms.string('bphWriteSpecificDecay:k0Fitted:bphAnalysis'),
    l0CandsLabel = cms.string('bphWriteSpecificDecay:l0Fitted:bphAnalysis'),
    b0CandsLabel = cms.string('bphWriteSpecificDecay:b0Fitted:bphAnalysis'),
    lbCandsLabel = cms.string('bphWriteSpecificDecay:lbFitted:bphAnalysis')
)

process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string('reco.root'),
    outputCommands = cms.untracked.vstring(
      "keep *_bphWriteSpecificDecay_*_*",
      "keep *_TriggerResults_*_HLT",
      "keep *_offlineSlimmedPrimaryVertices_*_*"
    ),
)

process.p = cms.Path(
    process.bphWriteSpecificDecay
    *
    process.bphHistoSpecificDecay
)

process.e = cms.EndPath(process.out)

