#include "AnalysisDataFormats/TauAnalysis/interface/NSVfitResonanceHypothesis.h"

NSVfitResonanceHypothesis::NSVfitResonanceHypothesis(const NSVfitResonanceHypothesis& bluePrint)
  : NSVfitResonanceHypothesisBase(bluePrint),
    eventHyp_(bluePrint.eventHyp_),
    p4_(bluePrint.p4_),
    dp4_(bluePrint.dp4_)
{
  size_t  numDaughters = daughters_.size();
  for ( size_t iDaughter = 0; iDaughter < numDaughters; iDaughter++ ) {
    this->daughter(iDaughter)->setMother(this);
  }
}

NSVfitResonanceHypothesis& NSVfitResonanceHypothesis::operator=(const NSVfitResonanceHypothesis& bluePrint)
{
  NSVfitResonanceHypothesisBase::operator=(bluePrint);
  eventHyp_ = bluePrint.eventHyp_;
  p4_ = bluePrint.p4_;
  dp4_ = bluePrint.dp4_;
  size_t  numDaughters = daughters_.size();
  for ( size_t iDaughter = 0; iDaughter < numDaughters; iDaughter++ ) {
    this->daughter(iDaughter)->setMother(this);
  }
  return (*this);
}

