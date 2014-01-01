#ifndef Bianchi_TutorialDAS2012_DiTauPairSelector_h
#define Bianchi_TutorialDAS2012_DiTauPairSelector_h

#include "DataFormats/Common/interface/RefVector.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "AnalysisDataFormats/TauAnalysis/interface/CompositePtrCandidateT1T2MEt.h"

#include <vector>

typedef SingleObjectSelector<
  std::vector<PATMuTauPair>,
  StringCutObjectSelector<PATMuTauPair> >  MuTauPairSelector;

typedef SingleObjectSelector<
  std::vector<PATMuTauPair>,
  StringCutObjectSelector<PATMuTauPair>,
  edm::RefVector<std::vector<PATMuTauPair>  > >  MuTauPairRefSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuTauPairSelector);
DEFINE_FWK_MODULE(MuTauPairRefSelector);

#endif
