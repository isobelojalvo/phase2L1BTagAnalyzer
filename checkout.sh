#!/bin/bash
scramv1 project CMSSW  CMSSW_10_1_5
cd CMSSW_10_1_5/src
eval `scramv1 runtime -sh`
git cms-init
git remote add cms-l1t-offline git@github.com:cms-l1t-offline/cmssw.git
git fetch cms-l1t-offline phase2-l1t-integration-CMSSW_10_1_5
git cms-merge-topic -u cms-l1t-offline:l1t-phase2-v2.16.8
#
# Tracklet Tracks
#
git remote add rekovic git@github.com:rekovic/cmssw.git
git fetch rekovic Tracklet-10_1_0_pre3
git cms-merge-topic -u rekovic:Tracklet-10_1_0_pre3-from-skinnari

git cms-addpkg L1Trigger/L1TCommon

git clone git@github.com:isobelojalvo/phase2L1BTagAnalyzer.git

scram b -j 8
