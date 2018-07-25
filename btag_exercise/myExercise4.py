#================================================
# Exercise 4 from CMSDAS 2018 B-tagging exercise
# 
# Instructions:
# Run this at the terminal before running the
# script, as it is necessary for accessing the 
# miniAOD files:
#
# % voms-proxy-init -voms cms
#
# Run the script: % python myExercise4.py
# 
# A .pdf output file will be generated.
#================================================

# Import necessary libraries/modules
import ROOT
import rootpy #hands down, a better version of PyROOT
import rootpy.plotting as plt
import pprint
from DataFormats.FWLite import Events, Handle
import pandas as pd
import numpy as np
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import matplotlib as mpl
import matplotlib.pyplot as pyplt
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.multiclass import OneVsRestClassifier
import sklearn

# Read input miniAOD files and convert the information we want into a 
# pandas.DataFrame

data = []
files = [
    #QCD
    '/store/relval/CMSSW_9_2_2/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_92X_upgrade2017_realistic_v1-v1/10000/14008288-C84D-E711-9EFD-0025905B85BC.root',
    #'/store/relval/CMSSW_9_2_2/RelValQCD_FlatPt_15_3000HS_13/MINIAODSIM/PU25ns_92X_upgrade2017_realistic_v1-v1/10000/30858183-C84D-E711-AB11-0CC47A7C3434.root',
    #TTbar
    '/store/relval/CMSSW_9_2_2/RelValTTbar_13/MINIAODSIM/PU25ns_92X_upgrade2017_realistic_v1-v1/10000/8E7EE25F-294E-E711-A5CC-0025905B8610.root',
    #'/store/relval/CMSSW_9_2_2/RelValTTbar_13/MINIAODSIM/PU25ns_92X_upgrade2017_realistic_v1-v1/10000/44D6F368-294E-E711-958E-0025905A612C.root',
]
events = Events(['root://xrootd-cms.infn.it/%s' % i for i in files])
handle = Handle('vector<pat::Jet>')
taggers = [
    ('JP', 'pfJetProbabilityBJetTags'),
    ('JPB', 'pfJetBProbabilityBJetTags'),
    ('SoftMu', 'softPFMuonBJetTags'),
    ('SoftEl', 'softPFElectronBJetTags'),
    ('CSV_IVF', 'pfCombinedInclusiveSecondaryVertexV2BJetTags'),
    ('CSV_AVR', 'pfCombinedSecondaryVertexV2BJetTags'),
    ('CvsL', 'pfCombinedCvsLJetTags'),
    ('CvsB', 'pfCombinedCvsBJetTags'),
    ('cMVA', 'pfCombinedMVAV2BJetTags'), #FOR REFERENCE ONLY
]
for event in events:
    event.getByLabel('slimmedJets', handle)
    jets = handle.product()
    for jet in jets:
        if jet.pt() < 20 or abs(jet.eta()) > 2.4: continue #basic selection
        #A more verbose, but more consistent version with dictionaries exists
        entry = [
            jet.pt(),
            jet.eta(),
            jet.mass(),
            abs(jet.hadronFlavour()), #5 - b-jet, 4 - c-jet, 0 - light
        ]
        entry.extend([jet.bDiscriminator(i) for _, i in taggers])
        data.append(entry)

# Convert the data and look at it
data = pd.DataFrame(
    data, 
    columns=(['pt', 'eta', 'mass', 'flavour'] + [i for i, _ in taggers])
)
#data.head()

# Check if there are columns with infs or NaNs, set to 0
for column in data.columns:
    if np.isinf(data[column]).any():
        print column, 'contains infs'
    if np.isnan(data[column]).any():
        print column, 'contains nans'

data.loc[np.isinf(data.SoftEl), 'SoftEl'] = 0
data.loc[np.isinf(data.SoftMu), 'SoftMu'] = 0

# We will now have to set the true labels for the training/testing of our discriminator. We will create two sets of labels:
# Binary labels: simply define if the jet is a b-jet or not
# multiclass labels: define three possible options, light, charm and b. It is the same information we have in the flavour column, but encoded differently

data['binary_target'] = 0
data.loc[(data.flavour == 5), 'binary_target'] = 1

data['isL'] = (data.flavour == 0)
data['isB'] = (data.flavour == 5)
data['isC'] = (data.flavour == 4)
data['clf_binary'] = 0
data['clf_multiclass'] = 0
data['clf_multiclass_C'] = 0

# As in every machine learning exercise, we need to split our data into two samples, a training and a testing.
train, test = train_test_split(data, test_size=0.3, random_state=42)

# We can now define our models. We will train two different BDTs, one trained on binary labels, one with multiclass labels.
features = [i for i, _ in taggers if i != 'cMVA'] + ['mass']
clf_binary = GradientBoostingClassifier(
    learning_rate=0.01, n_estimators=1000, 
    subsample=0.8, random_state=13,
    min_samples_leaf=int(0.01*len(train)),
    max_depth=5,
    verbose=1,
)
clf_multiclass = OneVsRestClassifier(sklearn.base.clone(clf_binary))

# Training, takes a while
clf_multiclass.fit(train[features], train[['isL', 'isC', 'isB']].as_matrix())
clf_binary.fit(train[features], train.binary_target)

# We can now store the output of the taggers in the columns we created before.
# Do not mind the warning coming from the command.
# You can see that the binary classifier outputs a N time 2 matrix corresponding to $p(non-b\, |\, jet)$, $p(b\, |\, jet)$, while the multiclass has three outputs $p(light\, |\, jet)$, $p(charm\, |\, jet)$, $p(b\, |\, jet)$.

test['clf_binary'] = clf_binary.predict_proba(test[features])[:,1]
test['clf_multiclass'] = clf_multiclass.predict_proba(test[features])[:,2]
test['clf_multiclass_C'] = clf_multiclass.predict_proba(test[features])[:,1]

# We can now print the ROC curves as in the previous exercise.
fig = pyplt.figure(figsize=(15, 15), dpi= 80, facecolor='w', edgecolor='k')
for algo, color in zip([
    'cMVA', 
    'clf_binary', 
    'clf_multiclass',
    ], 'rgb'):    
    for bkg, style in zip([4, 0], ['-', '--']):
        mask = (test.flavour != bkg)
        jets = test[mask]
        fakes_positive_rate, true_positive_rate, _ = roc_curve(jets.isB, jets[algo])
        pyplt.plot(true_positive_rate, fakes_positive_rate, '%s%s' % (color, style))

import matplotlib.lines as mlines
pyplt.legend(
    loc='best',
    handles=[
        mlines.Line2D([],[], color='red', ls='-', label='cMVAv2'),
        mlines.Line2D([],[], color='green', ls='-', label='clf_binary'),
        mlines.Line2D([],[], color='blue', ls='-', label='clf_multiclass'),
        mlines.Line2D([],[], color='k', ls='-', label='b vs. light'),
        mlines.Line2D([],[], color='k', ls='--', label='b vs. charm'),        
        ])
        
pyplt.xlabel('efficiency')
pyplt.ylabel('fake rate')
pyplt.grid(True)
pyplt.yscale('log', nonposy='clip')
#pyplt.show()

pyplt.savefig("exercise4_roc.pdf")
