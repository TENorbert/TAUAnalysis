#! /bin/sh

echo 'Creating cfg for DY->ll'
cp ./patTuple_PAT_MuTau_cfg.py ./patTuple_PAT_MuTau_DYJets_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g'                                               ./patTuple_PAT_MuTau_cfg.py   >> ./patTuple_PAT_MuTau_DYJets_cfg_1.py
sed 's/sample      = "DYJets"/sample      = "DYJets"/g'                                       ./patTuple_PAT_MuTau_DYJets_cfg_1.py >> ./patTuple_PAT_MuTau_DYJets_cfg_2.py
sed 's/process.source.fileNames = fileListDYJets/process.source.fileNames = fileListDYJets/g' ./patTuple_PAT_MuTau_DYJets_cfg_2.py >> ./patTuple_PAT_MuTau_DYJets_cfg_3.py
mv ./patTuple_PAT_MuTau_DYJets_cfg_3.py ./patTuple_PAT_MuTau_DYJets_cfg.py
rm ./patTuple_PAT_MuTau_DYJets_cfg_*.py

echo 'Creating cfg for W+jets'
cp ./patTuple_PAT_MuTau_cfg.py ./patTuple_PAT_MuTau_WJets_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g'                                               ./patTuple_PAT_MuTau_cfg.py   >> ./patTuple_PAT_MuTau_WJets_cfg_1.py
sed 's/sample      = "DYJets"/sample      = "WJets"/g'                                       ./patTuple_PAT_MuTau_WJets_cfg_1.py >> ./patTuple_PAT_MuTau_WJets_cfg_2.py
sed 's/process.source.fileNames = fileListDYJets/process.source.fileNames = fileListWJets/g' ./patTuple_PAT_MuTau_WJets_cfg_2.py >> ./patTuple_PAT_MuTau_WJets_cfg_3.py
mv ./patTuple_PAT_MuTau_WJets_cfg_3.py ./patTuple_PAT_MuTau_WJets_cfg.py
rm ./patTuple_PAT_MuTau_WJets_cfg_*.py

echo 'Creating cfg for ttbar'
cp ./patTuple_PAT_MuTau_cfg.py ./patTuple_PAT_MuTau_TTJets_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g'                                               ./patTuple_PAT_MuTau_cfg.py   >> ./patTuple_PAT_MuTau_TTJets_cfg_1.py
sed 's/sample      = "DYJets"/sample      = "TTJets"/g'                                       ./patTuple_PAT_MuTau_TTJets_cfg_1.py >> ./patTuple_PAT_MuTau_TTJets_cfg_2.py
sed 's/process.source.fileNames = fileListDYJets/process.source.fileNames = fileListTTJets/g' ./patTuple_PAT_MuTau_TTJets_cfg_2.py >> ./patTuple_PAT_MuTau_TTJets_cfg_3.py
mv ./patTuple_PAT_MuTau_TTJets_cfg_3.py ./patTuple_PAT_MuTau_TTJets_cfg.py
rm ./patTuple_PAT_MuTau_TTJets_cfg_*.py

echo 'Creating cfg for qq->H(130)->tautau'
cp ./patTuple_PAT_MuTau_cfg.py ./patTuple_PAT_MuTau_VBFH130_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g'                                               ./patTuple_PAT_MuTau_cfg.py   >> ./patTuple_PAT_MuTau_VBFH130_cfg_1.py
sed 's/sample      = "DYJets"/sample      = "VBFH130"/g'                                       ./patTuple_PAT_MuTau_VBFH130_cfg_1.py >> ./patTuple_PAT_MuTau_VBFH130_cfg_2.py
sed 's/process.source.fileNames = fileListDYJets/process.source.fileNames = fileListVBFH130/g' ./patTuple_PAT_MuTau_VBFH130_cfg_2.py >> ./patTuple_PAT_MuTau_VBFH130_cfg_3.py
mv ./patTuple_PAT_MuTau_VBFH130_cfg_3.py ./patTuple_PAT_MuTau_VBFH130_cfg.py
rm ./patTuple_PAT_MuTau_VBFH130_cfg_*.py

echo 'Creating cfg for gg->H(130)->tautau'
cp ./patTuple_PAT_MuTau_cfg.py ./patTuple_PAT_MuTau_GGFH130_cfg.py
sed 's/runOnMC     = True/runOnMC     = True/g'                                               ./patTuple_PAT_MuTau_cfg.py   >> ./patTuple_PAT_MuTau_GGFH130_cfg_1.py
sed 's/sample      = "DYJets"/sample      = "GGFH130"/g'                                       ./patTuple_PAT_MuTau_GGFH130_cfg_1.py >> ./patTuple_PAT_MuTau_GGFH130_cfg_2.py
sed 's/process.source.fileNames = fileListDYJets/process.source.fileNames = fileListGGFH130/g' ./patTuple_PAT_MuTau_GGFH130_cfg_2.py >> ./patTuple_PAT_MuTau_GGFH130_cfg_3.py
mv ./patTuple_PAT_MuTau_GGFH130_cfg_3.py ./patTuple_PAT_MuTau_GGFH130_cfg.py
rm ./patTuple_PAT_MuTau_GGFH130_cfg_*.py

echo 'Creating cfg for data'
cp ./patTuple_PAT_MuTau_cfg.py ./patTuple_PAT_MuTau_Data_cfg.py
sed 's/runOnMC     = True/runOnMC     = False/g'                                               ./patTuple_PAT_MuTau_cfg.py   >> ./patTuple_PAT_MuTau_Data_cfg_1.py
sed 's/sample      = "DYJets"/sample      = "Data"/g'                                       ./patTuple_PAT_MuTau_Data_cfg_1.py >> ./patTuple_PAT_MuTau_Data_cfg_2.py
sed 's/process.source.fileNames = fileListDYJets/process.source.fileNames = fileListData/g' ./patTuple_PAT_MuTau_Data_cfg_2.py >> ./patTuple_PAT_MuTau_Data_cfg_3.py
mv ./patTuple_PAT_MuTau_Data_cfg_3.py ./patTuple_PAT_MuTau_Data_cfg.py
rm ./patTuple_PAT_MuTau_Data_cfg_*.py
