import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask

process = cms.Process("bphAnalysis")

patAlgosToolsTask = getPatAlgosToolsTask(process)

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
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
    fileNames = cms.untracked.vstring('file:/home/ronchese/data/store/data/Run2017C/Charmonium/AOD/PromptReco-v2/000/300/079/00000/1C712689-2177-E711-944B-02163E01472F.root')
#    fileNames = cms.untracked.vstring('file:/home/ronchese/data/store/data/Run2017C/Charmonium/MINIAOD/PromptReco-v1/000/299/368/00000/64EE2B4F-866D-E711-A9DF-02163E01A6CE.root')
#    fileNames = cms.untracked.vstring('file:/lustre/cmswork/cmsdata/ronchese/store/data/Run2017B/Charmonium/MINIAOD/PromptReco-v2/000/298/853/00000/0C3845CC-8A68-E711-AAC9-02163E0118E2.root')
)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v2', '')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
patAlgosToolsTask.add(process.patCandidatesTask)
process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
patAlgosToolsTask.add(process.selectedPatCandidatesTask)
process.load('PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff')
patAlgosToolsTask.add(process.cleanPatCandidatesTask)

### add trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, None, None, None, None, '' )

process.selectedPatMuons.cut = cms.string('muonID(\"TMOneStationTight\")'
                    ' && abs(innerTrack.dxy) < 0.3'
                    ' && abs(innerTrack.dz)  < 20.'
                    ' && innerTrack.hitPattern.trackerLayersWithMeasurement > 5'
                    ' && innerTrack.hitPattern.pixelLayersWithMeasurement > 0'
                    ' && innerTrack.quality(\"highPurity\")'
                    )

#make patTracks
from PhysicsTools.PatAlgos.tools.trackTools import makeTrackCandidates
makeTrackCandidates(process,
                    label        = 'TrackCands',                  # output collection
                    tracks       = cms.InputTag('generalTracks'), # input track collection
                    particleType = 'pi+',                         # particle type (for assigning a mass)
                    preselection = 'pt > 0.7',                    # preselection cut on candidates
                    selection    = 'pt > 0.7',                    # selection on PAT Layer 1 objects
                    isolation    = {},                            # isolations to use (set to {} for None)
                    isoDeposits  = [],
                    mcAs         = None                           # replicate MC match as the one used for Muons
)
process.patTrackCands.embedTrack = True

# reduce MC genParticles a la miniAOD
process.load('PhysicsTools.PatAlgos.slimming.genParticles_cff')
process.packedGenParticles.inputVertices = cms.InputTag('offlinePrimaryVertices')

from PhysicsTools.PatAlgos.tools.coreTools import runOnData
runOnData( process, outputModules = [] )

from BPHAnalysis.SpecificDecay.recoSelectForWrite_cfi import recoSelect

process.bphWriteSpecificDecay = cms.EDProducer('BPHWriteSpecificDecay',
    pVertexLabel = cms.string('offlinePrimaryVertices'),
    pfCandsLabel = cms.string('particleFlow::RECO'),
    patMuonLabel = cms.string('selectedPatMuons'),
    k0CandsLabel = cms.string('generalV0Candidates:Kshort:RECO'),
    l0CandsLabel = cms.string('generalV0Candidates:Lambda:RECO'),
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
    bcName     = cms.string('bcFitted'),
    x3872Name  = cms.string('x3872Fitted'),
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
    lbCandsLabel = cms.string('bphWriteSpecificDecay:lbFitted:bphAnalysis'),
    bcCandsLabel = cms.string('bphWriteSpecificDecay:bcFitted:bphAnalysis'),
    x3872CandsLabel = cms.string('bphWriteSpecificDecay:x3872Fitted:bphAnalysis')
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
    process.bphHistoSpecificDecay,
    patAlgosToolsTask
)

process.e = cms.EndPath(process.out)

