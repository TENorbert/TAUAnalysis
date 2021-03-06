#include "TauAnalysis/CandidateTools/plugins/NSVfitTauToHadLikelihoodPhaseSpace.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitAlgorithmBase.h"

#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitTauToHadHypothesis.h"

#include <TMath.h>

NSVfitTauToHadLikelihoodPhaseSpace::NSVfitTauToHadLikelihoodPhaseSpace(const edm::ParameterSet& cfg)
  : NSVfitSingleParticleLikelihood(cfg)
{
// nothing to be done yet...
}

NSVfitTauToHadLikelihoodPhaseSpace::~NSVfitTauToHadLikelihoodPhaseSpace()
{
// nothing to be done yet...
}

void NSVfitTauToHadLikelihoodPhaseSpace::beginJob(NSVfitAlgorithmBase* algorithm)
{
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_visEnFracX, pluginName_);
  algorithm->requestFitParameter(prodParticleLabel_, nSVfit_namespace::kTau_phi_lab,    pluginName_);
}

double NSVfitTauToHadLikelihoodPhaseSpace::operator()(const NSVfitSingleParticleHypothesis* hypothesis) const
{
//--- compute negative log-likelihood for tau lepton decay "leg"
//    to be compatible with three-body decay,
//    assuming constant matrix element,
//    so that energy and angular distribution of decay products is solely determined by phase-space
//
//    NOTE: the parametrization of the three-body decay phase-space is taken from the PDG:
//          K. Nakamura et al. (Particle Data Group), J. Phys. G 37, 075021 (2010);
//          formula 38.20a
//
  //std::cout << "<NSVfitTauToHadLikelihoodPhaseSpace::operator()>:" << std::endl;

  const NSVfitTauToHadHypothesis* hypothesis_T = dynamic_cast<const NSVfitTauToHadHypothesis*>(hypothesis);
  assert(hypothesis_T != 0);

  if ( this->verbosity_ ) std::cout << "<NSVfitTauToHadLikelihoodPhaseSpace::operator()>:" << std::endl;

  double decayAngle = hypothesis_T->decay_angle_rf();

  if ( this->verbosity_ ) std::cout << " decayAngle = " << decayAngle << std::endl;

  double prob = TMath::Sin(decayAngle);
  
  double nll = 0.;
  if ( prob > 0. ) {
    nll = -TMath::Log(prob);
  } else {
    nll = std::numeric_limits<float>::max();
  }
  
  if ( this->verbosity_ ) std::cout << "--> nll = " << nll << std::endl;

  return nll;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_EDM_PLUGIN(NSVfitSingleParticleLikelihoodPluginFactory, NSVfitTauToHadLikelihoodPhaseSpace, "NSVfitTauToHadLikelihoodPhaseSpace");
