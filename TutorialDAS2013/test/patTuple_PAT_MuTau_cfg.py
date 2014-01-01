# import the process definition and the main PAT sequences
from PhysicsTools.PatAlgos.patTemplate_cfg import *

# load EDServices, Geometry and Jet Energy Corrections from the data base
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

# name to be appended to some PF-based collections
postfix     = "PFlow"
# a flag to enable/disable MC matching: DON'T CHANGE THE SPACING
runOnMC     = True
# name to be appended to the output file name: DON'T CHANGE THE SPACING
sample      = "DYJets"

from Configuration.PyReleaseValidation.autoCond import autoCond
process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )

if runOnMC:
    process.GlobalTag.globaltag = cms.string('START42_V13::All')

else:
    process.GlobalTag.globaltag = cms.string('GR_R_42_V19::All')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.fileNames = cms.untracked.vstring(
    '/store/user/bianchi/TutorialDAS2012/DYJets-MuTau-50-madgraph-PUS6-Pisa/patTuples_10_1_BUN.root'
    )


#################### input file names ###############
from Bianchi.TutorialDAS2012.DYJets_cff  import fileList as fileListDYJets
from Bianchi.TutorialDAS2012.WJets_cff   import fileList as fileListWJets
from Bianchi.TutorialDAS2012.TTJets_cff  import fileList as fileListTTJets
from Bianchi.TutorialDAS2012.VBFH130_cff import fileList as fileListVBFH130
from Bianchi.TutorialDAS2012.GGFH130_cff import fileList as fileListGGFH130
from Bianchi.TutorialDAS2012.Data_cff    import fileList as fileListData

### IF NEEDED COMMENT IT BUT DON'T CHANGE THE SPACING
process.source.fileNames = fileListDYJets

#process.source.eventsToProcess = cms.untracked.VEventRange(
#    '1:751063'
#    )


################### counter  ########################

process.allEventsFilter = cms.EDFilter(
    "AllEventsFilter"
    )

################### gen listing  ####################

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree1 = cms.EDAnalyzer(
    "ParticleListDrawer",
    src = cms.InputTag("genParticles"),
    maxEventsToPrint  = cms.untracked.int32(1)
    )


################### vertex sequence #################

process.selectedPrimaryVertices = cms.EDFilter(
    "VertexSelector",
    src = cms.InputTag('offlinePrimaryVertices'),
    cut = cms.string("isValid & ndof >= 4 & z > -24 & z < +24 & position.Rho < 2."),
    filter = cms.bool(False)                                          
)

process.primaryVertexCounter = cms.EDFilter(
    "VertexCountFilter",
    src = cms.InputTag('selectedPrimaryVertices'),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

################### jet sequence ####################

process.load('RecoJets.Configuration.RecoPFJets_cff')

process.kt6PFJets.doRhoFastjet  = True
process.kt6PFJets.Rho_EtaMax    = cms.double(4.4)
process.kt6PFJets.Ghost_EtaMax  = cms.double(5.0)
process.ak5PFJets.doAreaFastjet = True

## re-run kt4PFJets within lepton acceptance to compute rho
process.load('RecoJets.JetProducers.kt4PFJets_cfi')
process.fjSequence = cms.Sequence(process.kt6PFJets+process.ak5PFJets)


################### pat specific ####################

from PhysicsTools.PatAlgos.tools.helpers import massSearchReplaceAnyInputTag

from PhysicsTools.PatAlgos.tools.coreTools import *
if not runOnMC:
    removeMCMatching(process,["All"])
    
removeSpecificPATObjects(process, ['Photons'],
                         outputInProcess=False)
removeCleaning(process,
               outputInProcess=False)

restrictInputToAOD(process, ['All'])

# add particle-flow isolation to the muons
from Bianchi.TutorialDAS2012.customizePAT  import *
addSelectedPFlowParticle(process)

# add particle-flow MET
from PhysicsTools.PatAlgos.tools.metTools import *
addPfMET(process, postfix)

# switch to aKt particle-flow jets with R=0.5
from PhysicsTools.PatAlgos.tools.jetTools import *
switchJetCollection(process,cms.InputTag('ak5PFJets'),
                    doJTA        = True,
                    doBTagging   = True,
                    jetCorrLabel = ('AK5PF', ['L2Relative', 'L3Absolute',]),
                    doType1MET   = False,
                    genJetCollection=cms.InputTag("ak5GenJets"),
                    doJetID      = True,
                    jetIdLabel   = 'ak5'
                    )

# add jet energy corrections
JEClevels = cms.vstring(['L2Relative', 'L3Absolute'])
if runOnMC:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute']
else:
    JEClevels = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

process.patJetCorrFactors.levels = JEClevels
process.patJetCorrFactors.rho    = cms.InputTag('kt6PFJets','rho')
process.patJetCorrFactors.useRho = True
process.patJets.jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactors") )


#################### tau sequence #######################

process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")

# use the HPS tau algorithm
from PhysicsTools.PatAlgos.tools.tauTools import *
switchToPFTauHPS(process, 
                 pfTauLabelOld = 'shrinkingConePFTauProducer',
                 pfTauLabelNew = 'hpsPFTauProducer'
                 )

# embed these data into the dataformat
getattr(process,"patTaus").embedIsolationTracks             = cms.bool(True)
getattr(process,"patTaus").embedSignalTracks                = cms.bool(True)
getattr(process,"patTaus").embedGenMatch                    = cms.bool(True)
getattr(process,"patTaus").embedLeadTrack                   = cms.bool(True)
getattr(process,"patTaus").embedLeadPFCand                  = cms.bool(True)
getattr(process,"patTaus").embedLeadPFChargedHadrCand       = cms.bool(True)
getattr(process,"patTaus").embedLeadPFNeutralCand           = cms.bool(True)
getattr(process,"patTaus").embedSignalPFCands               = cms.bool(True)
getattr(process,"patTaus").embedSignalPFChargedHadrCands    = cms.bool(True)
getattr(process,"patTaus").embedSignalPFNeutralHadrCands    = cms.bool(True)
getattr(process,"patTaus").embedSignalPFGammaCands          = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFCands            = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFChargedHadrCands = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFNeutralHadrCands = cms.bool(True)
getattr(process,"patTaus").embedIsolationPFGammaCands       = cms.bool(True)
getattr(process,"patTaus").embedGenJetMatch                 = cms.bool(True)

# configure reco-tau-gen matching
process.tauMatch.maxDeltaR                = 0.15
process.tauMatch.resolveAmbiguities       = cms.bool(False)
process.tauGenJetMatch.resolveAmbiguities = cms.bool(False)
process.tauGenJetMatch.maxDeltaR          = 0.15
process.tauGenJetMatch.maxDPtRel          = 999


########################  pat::muon  #############################

# defined in customize.py
addPFMuonIsolation(process,process.patMuons)
getattr(process,"patMuons").embedTrack = True

######################## pat::electron ###########################

# defined in customize.py
addPFElectronIsolation(process,process.patElectrons)
getattr(process,"patElectrons").embedTrack    = True
getattr(process,"patElectrons").embedGsfTrack = True

######################## pat::jet ################################

# discard jets below 10 GeV
getattr(process,"selectedPatJets").cut = cms.string('pt>10 && abs(eta)<5.0')

######################## pat::trigger ############################

from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process )
process.patTriggerEvent.processName = '*'

if hasattr(process,"patTrigger"):
    process.patTrigger.processName = '*'
    

######################## embedding ###############################

# create a new collection of patMuons with new informations
process.selectedPatMuonsUserEmbedded = cms.EDProducer(
    "MuonsUserEmbedded",
    muonTag   = cms.InputTag("selectedPatMuons"),
    vertexTag = cms.InputTag("offlinePrimaryVertices")
    )
# create a new collection of patTaus with new informations
process.selectedPatTausUserEmbedded = cms.EDProducer(
    "TausUserEmbedded",
    tauTag    = cms.InputTag("selectedPatTaus"),
    vertexTag = cms.InputTag("offlinePrimaryVertices"),
    )

####################### muon selection ###########################

process.muPtEtaID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.1"+
                     " && isTrackerMuon && isGlobalMuon"+
                     " && numberOfMatches>=2"+
                     " && globalTrack.isNonnull "+
                     " && globalTrack.hitPattern.numberOfValidMuonHits>=1"+
                     " && globalTrack.hitPattern.numberOfValidPixelHits>=1"+
                     " && globalTrack.hitPattern.numberOfValidTrackerHits>=10"+
                     " && globalTrack.normalizedChi2<10"+
                     " && globalTrack.ptError/globalTrack.pt<0.1"+
                     " && abs(userFloat('dxyWrtPV'))<0.045 && abs(userFloat('dzWrtPV'))<0.2"
                     ),
    filter = cms.bool(False)
    )

# collection of loosely identified muons
process.muPtEtaRelID = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("selectedPatMuonsUserEmbedded"),
    cut = cms.string("pt>15 && abs(eta)<2.4 && isGlobalMuon"),
    filter = cms.bool(False)
    )

####################### tau selection #############################

process.tauPtEtaIDAgMuAgElec  = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("selectedPatTausUserEmbedded"),
    cut = cms.string("pt>19 && abs(eta)<2.3"+
                     " && tauID('decayModeFinding')>0.5"+
                     " && tauID('byLooseCombinedIsolationDeltaBetaCorr')>0.5"+
                     " && userFloat('dzWrtPV')<0.2"+
                     " && tauID('againstMuonTight')>0.5"+
                     " && tauID('againstElectronLoose')>0.5"),
    filter = cms.bool(False)
    )

####################### pairing ##################################

# require the muon an ta to be separated by more than dR=0.5
process.atLeastOneMuTau = cms.EDProducer(
    "CandViewShallowCloneCombiner",
    decay = cms.string("muPtEtaID tauPtEtaIDAgMuAgElec"),
    cut = cms.string("sqrt((daughter(0).eta-daughter(1).eta)*(daughter(0).eta-daughter(1).eta)+  min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi)  ) *  min( abs(daughter(0).phi-daughter(1).phi), 2*3.1415926 - abs(daughter(0).phi-daughter(1).phi)  )  )>0.5"),
    checkCharge = cms.bool(False)
    )
process.atLeastOneMuTauCounter = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("atLeastOneMuTau"),
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999),
    )

###################### final sequences ##############################

process.atLeastOneGoodVertexSequence = cms.Sequence(
    process.selectedPrimaryVertices*
    process.primaryVertexCounter
    )

process.alLeastOneMuTauSequence = cms.Sequence(
    process.atLeastOneMuTau*
    process.atLeastOneMuTauCounter
    )

process.muLegSequence = cms.Sequence(
    process.muPtEtaID+process.muPtEtaRelID
    )

process.tauLegSequence = cms.Sequence(
    process.tauPtEtaIDAgMuAgElec
    )

########################## di-tau  ##########################

# configure the di-tau sequence <=> create a collection of di-tau objects with
# tau-pair mass informations
process.load("Bianchi.TutorialDAS2012.diTauReconstruction_cff")
process.diTau = process.allMuTauPairs.clone()
process.diTau.srcLeg1  = cms.InputTag("muPtEtaID")
process.diTau.srcLeg2  = cms.InputTag("tauPtEtaIDAgMuAgElec")
process.diTau.srcMET   = cms.InputTag("patMETsPFlow")
process.diTau.dRmin12  = cms.double(0.5)
process.diTau.doSVreco = cms.bool(True)

if not runOnMC:
    process.diTau.srcGenParticles = ""

process.diTauSequence = cms.Sequence( process.diTau )


########################## analyzer  ##########################
# refer to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolPhysicsWithTausExercise
# for the exact meaning

process.muTauAnalyzer = cms.EDAnalyzer(
    "MuTauAnalyzer",
    diTaus         = cms.InputTag("diTau"),
    jets           = cms.InputTag("selectedPatJets"),
    rawMet         = cms.InputTag("patMETsPFlow"),
    muons          = cms.InputTag("muPtEtaID"),
    muonsRel       = cms.InputTag("muPtEtaRelID"),
    vertices       = cms.InputTag("selectedPrimaryVertices"),
    triggerResults = cms.InputTag("patTriggerEvent"),
    isMC           = cms.bool(runOnMC),
    deltaRLegJet   = cms.untracked.double(0.5),
    minCorrPt      = cms.untracked.double(15.),
    minJetID       = cms.untracked.double(0.5),
    verbose        = cms.untracked.bool( False ),
    )

########################## path ###############################
# refer to https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolPhysicsWithTausExercise
# for the exact meaning

process.skim = cms.Sequence(
    process.allEventsFilter+
    process.atLeastOneGoodVertexSequence*
    process.PFTau*
    process.fjSequence*
    process.patDefaultSequence*
    process.selectedPatMuonsUserEmbedded*
    process.selectedPatTausUserEmbedded*
    process.muLegSequence*
    process.tauLegSequence*
    process.alLeastOneMuTauSequence*
    process.diTauSequence*
    process.muTauAnalyzer*
    process.printTree1
    )

# use the filtered "selectedPrimaryVertices" vertices instead of the
# default collection
massSearchReplaceAnyInputTag(process.skim,
                             "offlinePrimaryVertices",
                             "selectedPrimaryVertices",
                             verbose=False)
process.selectedPrimaryVertices.src = cms.InputTag('offlinePrimaryVertices')


if not runOnMC:
    process.skim.remove(process.printTree1)

# the final PATH!
process.p = cms.Path(process.skim)


########################## output ###############################

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("treeMuTau_"+sample+".root")
    )

from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out.outputCommands = cms.untracked.vstring('drop *',
                                                   *patEventContentNoCleaning ) 
from Bianchi.TutorialDAS2012.eventContent_cff import eventContent
process.out.outputCommands.extend( eventContent )

process.out.SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('p')
    )
process.out.fileName = cms.untracked.string('patTuples_MuTau_'+sample+'.root')
process.outpath = cms.EndPath()

