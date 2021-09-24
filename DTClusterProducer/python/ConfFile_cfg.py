import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi");
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi");
process.load("Geometry.CaloEventSetup.CaloTopology_cfi");
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = '120X_mcRun3_2021_realistic_v2'

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/home/users/mcitron/gen-sim-reco-test.root'
    )
)

process.demo = cms.EDProducer('DTClusterProducer',
)
process.output = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "~/tempOutput.root" ),
    )

process.p = cms.Path(process.demo)
process.Output = cms.EndPath(process.output)
