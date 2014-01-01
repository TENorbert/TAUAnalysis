import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.tools.jetTools import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.pfTools import *
from CommonTools.ParticleFlow.Isolation.tools_cfi import *

###################a#################################################
def addSelectedPFlowParticle(process,verbose=False):
    if verbose:
        print "[Info] Adding pf-particles (for pf-isolation and pf-seed pat-leptons)"
    process.load("CommonTools.ParticleFlow.ParticleSelectors.pfSortByType_cff")
    process.load("CommonTools.ParticleFlow.pfNoPileUp_cff")

    process.pfAllChargedCandidates = process.pfAllChargedHadrons.clone()
    process.pfAllChargedCandidates.pdgId.extend(process.pfAllMuons.pdgId.value())
    process.pfAllChargedCandidates.pdgId.extend(process.pfAllElectrons.pdgId.value())
    
    process.pfPU = cms.EDProducer(
        "TPPFCandidatesOnPFCandidates",
        enable =  cms.bool( True ),
        verbose = cms.untracked.bool( False ),
        name = cms.untracked.string("puPFCandidates"),
        topCollection = cms.InputTag("pfNoPileUp"),
        bottomCollection = cms.InputTag("particleFlow"),
        )
    process.pfAllChargedHadronsPU = process.pfAllChargedHadrons.clone(src='pfPU')

    #################################################################
    process.pfNoCharged = process.pfPU.clone(
        name = cms.untracked.string("noChargedPFCandidates"),
        topCollection = cms.InputTag("pfAllChargedHadrons"),
        bottomCollection = cms.InputTag("particleFlow"),
        )
    process.pfAllNeutral = cms.EDFilter(
        "PdgIdPFCandidateSelector",
        pdgId = cms.vint32(111, 130, 310, 2112, 22),
        src = cms.InputTag("pfNoCharged")
        )
    #################################################################

    process.pfCandidateSelectionByType = cms.Sequence(
        process.pfNoPileUpSequence *
        ( process.pfAllNeutralHadrons +
          process.pfAllChargedHadrons +
          process.pfAllChargedCandidates +
          process.pfAllPhotons +
          (process.pfPU * process.pfAllChargedHadronsPU )
          )  +
        process.pfAllMuons +
        process.pfAllElectrons +
        ( process.pfNoCharged+process.pfAllNeutral)
        )
    process.pfPileUp.Enable              = True
    process.pfPileUp.checkClosestZVertex = True
    process.pfPileUp.Vertices            = "offlinePrimaryVertices"
    process.pfAllMuons.src               = "pfNoPileUp"
    process.pfAllElectrons.src           = "pfNoPileUp"
    
    
###################a#################################################
def addPFMuonIsolation(process,module,postfix="",verbose=False):
    if verbose:
        print "[Info] Adding particle isolation to muon with postfix '"+postfix+"'"

    #if not hasattr(process, "pfCandidateSelectionByType"):
    #    addSelectedPFlowParticle(process,verbose=verbose)
        
    #setup correct src of isolated object
    setattr(process,"isoDepMuonWithCharged"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllChargedCandidates'))
    setattr(process,"isoDepMuonWithChargedHadrons"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllChargedHadrons'))
    setattr(process,"isoDepMuonWithMuons"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllMuons'))
    setattr(process,"isoDepMuonWithElectrons"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllElectrons'))
    setattr(process,"isoDepMuonWithChargedPU"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllChargedHadronsPU'))
    setattr(process,"isoDepMuonWithNeutral"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllNeutralHadrons'))
    setattr(process,"isoDepMuonWithPhotons"+postfix,
            isoDepositReplace(module.muonSource,
                              'pfAllPhotons'))

    #compute isolation values form deposits
    process.load("CommonTools.ParticleFlow.Isolation.pfMuonIsolationFromDeposits_cff")
    setattr(process,"isoValMuonWithChargedPU",
            process.isoValMuonWithCharged.clone())
    getattr(process,"isoValMuonWithChargedPU").deposits[0].src="isoDepMuonWithChargedPU"
    setattr(process,"isoValMuonWithChargedHadrons",
            process.isoValMuonWithCharged.clone())
    getattr(process,"isoValMuonWithChargedHadrons").deposits[0].src="isoDepMuonWithChargedHadrons"
    setattr(process,"isoValMuonWithMuons",
            process.isoValMuonWithCharged.clone())
    getattr(process,"isoValMuonWithMuons").deposits[0].src="isoDepMuonWithMuons"
    setattr(process,"isoValMuonWithElectrons",
            process.isoValMuonWithCharged.clone())
    getattr(process,"isoValMuonWithElectrons").deposits[0].src="isoDepMuonWithElectrons"

    if postfix!="":
        setattr(process,"isoValMuonWithCharged"+postfix,
                process.isoValMuonWithCharged.clone())
        getattr(process,"isoValMuonWithCharged"+postfix).deposits[0].src="isoDepMuonWithCharged"+postfix
        getattr(process,"isoValMuonWithMuons").deposits[0].src="isoDepMuonWithMuons"
        setattr(process,"isoValMuonWithElectrons",
                process.isoValMuonWithCharged.clone())
        getattr(process,"isoValMuonWithElectrons").deposits[0].src="isoDepMuonWithElectrons"
        setattr(process,"isoValMuonWithChargedHadrons"+postfix,
                process.isoValMuonWithCharged.clone())
        getattr(process,"isoValMuonWithChargedHadrons"+postfix).deposits[0].src="isoDepMuonWithChargedHadrons"+postfix
        setattr(process,"isoValMuonWithChargedPU"+postfix,
                process.isoValMuonWithChargedPU.clone())
        getattr(process,"isoValMuonWithChargedPU"+postfix).deposits[0].src="isoDepMuonWithChargedPU"+postfix
        setattr(process,"isoValMuonWithNeutral"+postfix,
                process.isoValMuonWithNeutral.clone())
        getattr(process,"isoValMuonWithNeutral"+postfix).deposits[0].src="isoDepMuonWithNeutral"+postfix
        setattr(process,"isoValMuonWithPhotons"+postfix,
                process.isoValMuonWithPhotons.clone())
        getattr(process,"isoValMuonWithPhotons"+postfix).deposits[0].src="isoDepMuonWithPhotons"+postfix
        
    setattr(process,"patMuonIsolationFromDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoValMuonWithCharged"+postfix) +
                         getattr(process,"isoValMuonWithChargedHadrons"+postfix) +
                         getattr(process,"isoValMuonWithMuons"+postfix) +
                         getattr(process,"isoValMuonWithElectrons"+postfix) +
                         getattr(process,"isoValMuonWithChargedPU"+postfix) +
                         getattr(process,"isoValMuonWithNeutral"+postfix) +
                         getattr(process,"isoValMuonWithPhotons"+postfix)
                         )
            )

    setattr(process,"patMuonIsoDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoDepMuonWithCharged"+postfix) +
                         getattr(process,"isoDepMuonWithChargedHadrons"+postfix) +
                         getattr(process,"isoDepMuonWithMuons"+postfix) +
                         getattr(process,"isoDepMuonWithElectrons"+postfix) +
                         getattr(process,"isoDepMuonWithChargedPU"+postfix) +
                         getattr(process,"isoDepMuonWithNeutral"+postfix) +
                         getattr(process,"isoDepMuonWithPhotons"+postfix)
                         )
            )
    setattr(process,"patMuonIsolationSequence"+postfix,
            cms.Sequence(getattr(process,"patMuonIsoDepositsSequence"+postfix) +
                         getattr(process,"patMuonIsolationFromDepositsSequence"+postfix)
                         )
            )
    
    getattr(process,"isoDepMuonWithCharged"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedCandidates"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithChargedHadrons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithMuons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllMuons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithElectrons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllElectrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithNeutral"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepMuonWithPhotons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )   
    getattr(process,"isoValMuonWithPhotons"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepMuonWithPhotons"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValMuonWithNeutral"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepMuonWithNeutral"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValMuonWithCharged"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepMuonWithCharged"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValMuonWithChargedPU"+postfix).deposits[0].vetos = getattr(
        process,"isoValMuonWithCharged"+postfix).deposits[0].vetos

    getattr(process,"isoValMuonWithChargedHadrons"+postfix).deposits[0].vetos = getattr(
        process,"isoValMuonWithCharged"+postfix).deposits[0].vetos

    getattr(process,"isoValMuonWithMuons"+postfix).deposits[0].vetos = getattr(
        process,"isoValMuonWithCharged"+postfix).deposits[0].vetos

    getattr(process,"isoValMuonWithElectrons"+postfix).deposits[0].vetos = getattr(
        process,"isoValMuonWithCharged"+postfix).deposits[0].vetos

    
    module.isoDeposits = cms.PSet(
        pfAllParticles   = cms.InputTag("isoDepMuonWithChargedPU"+postfix),
        pfChargedHadrons = cms.InputTag("isoDepMuonWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoDepMuonWithNeutral"+postfix),
        pfPhotons        = cms.InputTag("isoDepMuonWithPhotons"+postfix),
        user = cms.VInputTag(
        cms.InputTag("isoDepMuonWithChargedHadrons"+postfix),
        cms.InputTag("isoDepMuonWithMuons"+postfix),
        cms.InputTag("isoDepMuonWithElectrons"+postfix)
        )
        )
    module.isolationValues = cms.PSet(
        pfAllParticles   = cms.InputTag("isoValMuonWithChargedPU"+postfix),
        pfChargedHadrons = cms.InputTag("isoValMuonWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoValMuonWithNeutral"+postfix),
        pfPhotons        = cms.InputTag("isoValMuonWithPhotons"+postfix),
        user = cms.VInputTag(
        cms.InputTag("isoValMuonWithChargedHadrons"+postfix),
        cms.InputTag("isoValMuonWithMuons"+postfix),
        cms.InputTag("isoValMuonWithElectrons"+postfix)
        )
        )
    
    process.patDefaultSequence.replace(module,
                                       getattr(process,"patMuonIsolationSequence"+postfix)+
                                       module
                                       )
    process.patDefaultSequence.replace(process.isoDepMuonWithCharged,
                                       process.pfCandidateSelectionByType+
                                       process.isoDepMuonWithCharged)


#################################################################################
def addPFElectronIsolation(process,module,postfix="",verbose=False):
    if verbose:
        print "[Info] Adding particle isolation to electron with postfix '"+postfix+"'"

    #if not hasattr(process, "pfCandidateSelectionByType"):
    #    addSelectedPFlowParticle(process,verbose=verbose)
        
    #setup correct src of isolated object
    setattr(process,"isoDepElectronWithCharged"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllChargedCandidates'))
    setattr(process,"isoDepElectronWithChargedHadrons"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllChargedHadrons'))
    setattr(process,"isoDepElectronWithMuons"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllElectrons'))
    setattr(process,"isoDepElectronWithElectrons"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllElectrons'))
    setattr(process,"isoDepElectronWithChargedPU"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllChargedHadronsPU'))
    setattr(process,"isoDepElectronWithNeutral"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllNeutralHadrons'))
    setattr(process,"isoDepElectronWithPhotons"+postfix,
            isoDepositReplace(module.electronSource,
                              'pfAllPhotons'))

    #compute isolation values form deposits
    process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolationFromDeposits_cff")
    setattr(process,"isoValElectronWithChargedPU",
            process.isoValElectronWithCharged.clone())
    getattr(process,"isoValElectronWithChargedPU").deposits[0].src="isoDepElectronWithChargedPU"
    setattr(process,"isoValElectronWithChargedHadrons",
            process.isoValElectronWithCharged.clone())
    getattr(process,"isoValElectronWithChargedHadrons").deposits[0].src="isoDepElectronWithChargedHadrons"
    setattr(process,"isoValElectronWithMuons",
            process.isoValElectronWithCharged.clone())
    getattr(process,"isoValElectronWithMuons").deposits[0].src="isoDepElectronWithMuons"
    setattr(process,"isoValElectronWithElectrons",
            process.isoValElectronWithCharged.clone())
    getattr(process,"isoValElectronWithElectrons").deposits[0].src="isoDepElectronWithElectrons"


    #compute isolation values form deposits
    process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolationFromDeposits_cff")
    if postfix!="":
        setattr(process,"isoValElectronWithCharged"+postfix,
                process.isoValElectronWithCharged.clone())
        getattr(process,"isoValElectronWithCharged"+postfix).deposits[0].src="isoDepElectronWithCharged"+postfix
        setattr(process,"isoValElectronWithCharged"+postfix,
                process.isoValElectronWithCharged.clone())
        getattr(process,"isoValElectronWithCharged"+postfix).deposits[0].src="isoDepElectronWithCharged"+postfix
        getattr(process,"isoValElectronWithMuons").deposits[0].src="isoDepElectronWithMuons"
        setattr(process,"isoValElectronWithElectrons",
                process.isoValElectronWithCharged.clone())
        getattr(process,"isoValElectronWithElectrons").deposits[0].src="isoDepElectronWithElectrons"
        setattr(process,"isoValElectronWithChargedHadrons"+postfix,
                process.isoValElectronWithCharged.clone())
        getattr(process,"isoValElectronWithChargedHadrons"+postfix).deposits[0].src="isoDepElectronWithChargedHadrons"+postfix
        setattr(process,"isoValElectronWithChargedPU"+postfix,
                process.isoValElectronWithChargedPU.clone())
        getattr(process,"isoValElectronWithChargedPU"+postfix).deposits[0].src="isoDepElectronWithChargedPU"+postfix,
     
        setattr(process,"isoValElectronWithNeutral"+postfix,
                process.isoValElectronWithNeutral.clone())
        getattr(process,"isoValElectronWithNeutral"+postfix).deposits[0].src="isoDepElectronWithNeutral"+postfix
        setattr(process,"isoValElectronWithPhotons"+postfix,
                process.isoValElectronWithPhotons.clone())
        getattr(process,"isoValElectronWithPhotons"+postfix).deposits[0].src="isoDepElectronWithPhotons"+postfix
        
    setattr(process,"patElectronIsolationFromDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoValElectronWithCharged"+postfix) +
                         getattr(process,"isoValElectronWithChargedHadrons"+postfix) +
                         getattr(process,"isoValElectronWithMuons"+postfix) +
                         getattr(process,"isoValElectronWithElectrons"+postfix) +
                         getattr(process,"isoValElectronWithChargedPU"+postfix) +
                         getattr(process,"isoValElectronWithNeutral"+postfix) +
                         getattr(process,"isoValElectronWithPhotons"+postfix)
                         )
            )

    setattr(process,"patElectronIsoDepositsSequence"+postfix,
            cms.Sequence(getattr(process,"isoDepElectronWithCharged"+postfix) +
                         getattr(process,"isoDepElectronWithChargedHadrons"+postfix) +
                         getattr(process,"isoDepElectronWithMuons"+postfix) +
                         getattr(process,"isoDepElectronWithElectrons"+postfix) +
                         getattr(process,"isoDepElectronWithChargedPU"+postfix) +
                         getattr(process,"isoDepElectronWithNeutral"+postfix) +
                         getattr(process,"isoDepElectronWithPhotons"+postfix)
                         )
            )
    setattr(process,"patElectronIsolationSequence"+postfix,
            cms.Sequence(getattr(process,"patElectronIsoDepositsSequence"+postfix) +
                         getattr(process,"patElectronIsolationFromDepositsSequence"+postfix)
                         )
            )

    getattr(process,"isoDepElectronWithCharged"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedCandidates"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithChargedHadrons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllChargedHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithMuons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllMuons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithElectrons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllElectrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithNeutral"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllNeutralHadrons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    getattr(process,"isoDepElectronWithPhotons"+postfix).ExtractorPSet = cms.PSet(
        Diff_z = cms.double(99999.99),
        ComponentName = cms.string('CandViewExtractor'),
        DR_Max = cms.double(1.0),
        Diff_r = cms.double(99999.99),
        inputCandView = cms.InputTag("pfAllPhotons"),
        DR_Veto = cms.double(1e-05),
        DepositLabel = cms.untracked.string('')
        )
    
    getattr(process,"isoValElectronWithPhotons"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepElectronWithPhotons"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValElectronWithNeutral"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepElectronWithNeutral"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.5)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValElectronWithCharged"+postfix).deposits = cms.VPSet(
        cms.PSet(
        src = cms.InputTag("isoDepElectronWithCharged"+postfix),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('Threshold(0.0)'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
        )
        )
    getattr(process,"isoValElectronWithChargedPU"+postfix).deposits[0].vetos = getattr(
        process,"isoValElectronWithCharged"+postfix).deposits[0].vetos

    getattr(process,"isoValElectronWithChargedHadrons"+postfix).deposits[0].vetos = getattr(
        process,"isoValElectronWithCharged"+postfix).deposits[0].vetos

    getattr(process,"isoValElectronWithMuons"+postfix).deposits[0].vetos = getattr(
        process,"isoValElectronWithCharged"+postfix).deposits[0].vetos

    getattr(process,"isoValElectronWithElectrons"+postfix).deposits[0].vetos = getattr(
        process,"isoValElectronWithCharged"+postfix).deposits[0].vetos


    module.isoDeposits = cms.PSet(
        pfAllParticles = cms.InputTag("isoDepElectronWithChargedPU"+postfix),
        pfChargedHadrons = cms.InputTag("isoDepElectronWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoDepElectronWithNeutral"+postfix),
        pfPhotons = cms.InputTag("isoDepElectronWithPhotons"+postfix),
        user = cms.VInputTag(
        cms.InputTag("isoDepElectronWithChargedHadrons"+postfix),
        cms.InputTag("isoDepElectronWithMuons"+postfix),
        cms.InputTag("isoDepElectronWithElectrons"+postfix)
        )
        )
    module.isolationValues = cms.PSet(
        pfAllParticles = cms.InputTag("isoValElectronWithChargedPU"+postfix),
        pfChargedHadrons = cms.InputTag("isoValElectronWithCharged"+postfix),
        pfNeutralHadrons = cms.InputTag("isoValElectronWithNeutral"+postfix),
        pfPhotons = cms.InputTag("isoValElectronWithPhotons"+postfix),
        user = cms.VInputTag(
        cms.InputTag("isoValElectronWithChargedHadrons"+postfix),
        cms.InputTag("isoValElectronWithMuons"+postfix),
        cms.InputTag("isoValElectronWithElectrons"+postfix)
        )
        )
    
    process.patDefaultSequence.replace(module,
                                       getattr(process,"patElectronIsolationSequence"+postfix)+
                                       module
                                       )

    process.patDefaultSequence.replace(process.isoDepElectronWithCharged,
                                       process.pfCandidateSelectionByType+
                                       process.isoDepElectronWithCharged)

####################################################################


