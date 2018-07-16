import FWCore.ParameterSet.Config as cms
import string

process = cms.Process('RECODQM')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework
process.load("DQM.Integration.config.environment_cfi")
process.dqmEnv.subSystemFolder = "CTPPS"
process.dqmEnv.eventInfoFolder = "EventInfo"
process.dqmSaver.path = ""
process.dqmSaver.tag = "CTPPS"

# process.source = cms.Source("NewEventStreamFileReader",
#     fileNames = cms.untracked.vstring(
#     *(
#     'file:///afs/cern.ch/user/n/nminafra/Work/public/SampicCMS/run319190_ls0130_streamPhysicsTOTEM3_StorageManager.dat',
#     )
#     )
# )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    *(
    'root:///afs/cern.ch/user/n/nminafra/Work/public/SampicCMS/319190_testAOD.root',
    # '/store/data/Run2018B/TOTEM40/RECO/PromptReco-v2/000/319/190/00000/00A9F18A-B481-E811-9ECD-02163E014BBF.root ',
    )
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_HLT_v7', '')

# raw-to-digi conversion
process.load("EventFilter.CTPPSRawToDigi.ctppsRawToDigi_cff")

# local RP reconstruction chain with standard settings
process.load("RecoCTPPS.Configuration.recoCTPPS_cff")

# CTPPS DQM modules
process.load("DQM.CTPPS.ctppsDQM_cff")

# rechits production
process.load('RecoCTPPS.TotemRPLocal.totemTimingRecHits_cfi')
process.totemTimingRecHits.baselinePoints = cms.int32(6);
process.totemTimingRecHits.saturationLimit = cms.double(0.85);
process.totemTimingRecHits.cfdFraction = cms.double(0.3);
process.totemTimingRecHits.hysteresis = cms.double(10e-3);
process.totemTimingRecHits.smoothingPoints = cms.int32(6);
process.totemTimingRecHits.lowPassFrequency = cms.double(0.6);


process.simpleAnalyzer = cms.EDAnalyzer('SimpleAnalyzer',
 tagDigi = cms.InputTag("totemTimingRawToDigi", "TotemTiming"),
 tagRecHit = cms.InputTag("totemTimingRecHits"),
)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string('output.root')
)

process.path = cms.Path(
    process.simpleAnalyzer
)

process.schedule = cms.Schedule(
process.path
)
