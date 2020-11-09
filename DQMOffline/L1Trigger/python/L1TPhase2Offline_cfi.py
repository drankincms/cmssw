import FWCore.ParameterSet.Config as cms

from DQMServices.Core.DQMEDAnalyzer import DQMEDAnalyzer
l1tPhase2CorrelatorOfflineDQM = DQMEDAnalyzer(
    "L1TPhase2CorrelatorOffline",
    verbose   = cms.untracked.bool(False),

    phase2L1PfInputTag = cms.untracked.InputTag("l1pfCandidates:PF"),
    phase2L1PuppiInputTag = cms.untracked.InputTag("l1pfCandidates:Puppi"),
    genJetsInputTag = cms.untracked.InputTag("ak4GenJetsNoNu"),
    genParticlesInputTag = cms.untracked.InputTag("genParticles"),
    isParticleGun = cms.bool(False),
    objects = cms.PSet(
        L1PF = cms.VInputTag("l1pfCandidates:PF",),
        L1Puppi = cms.VInputTag("l1pfCandidates:Puppi",),
    ),

    histFolder = cms.string('L1T/L1TObjects/L1TPhase2/'),

    histDefinitions=cms.PSet(
        resVsPt=cms.PSet(
            name=cms.untracked.string("resVsPt"),
            title=cms.untracked.string("resVsPt"),
            nbinsX=cms.untracked.uint32(20),
            xmin=cms.untracked.double(0.),
            xmax=cms.untracked.double(100.),
        ),
        resVsEta=cms.PSet(
            name=cms.untracked.string("resVsEta"),
            title=cms.untracked.string("resVsEta"),
            nbinsX=cms.untracked.uint32(20),
            xmin=cms.untracked.double(-5.),
            xmax=cms.untracked.double(5.),
        ),
        ptDist=cms.PSet(
            name=cms.untracked.string("ptDist"),
            title=cms.untracked.string("ptDist"),
            nbinsX=cms.untracked.uint32(20),
            xmin=cms.untracked.double(0.),
            xmax=cms.untracked.double(100.),
        ),
        etaDist=cms.PSet(
            name=cms.untracked.string("etaDist"),
            title=cms.untracked.string("etaDist"),
            nbinsX=cms.untracked.uint32(20),
            xmin=cms.untracked.double(-5.),
            xmax=cms.untracked.double(5.),
        ),
    ),

)

l1tPhase2OfflineDQM = cms.Sequence(
                          l1tPhase2CorrelatorOfflineDQM
                          )
