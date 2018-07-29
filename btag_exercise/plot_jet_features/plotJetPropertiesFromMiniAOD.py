####################################################################
# Plot jet pT, eta, phi, mass, b-jet discriminators from a miniAOD
# file of a simulated t-tbar event, separating the dataset into 
# hadronFlavor == 5 (b-jets) and hadronFlavor != 5 (non-b).
# 
# Instructions:
# 1. Before running this script, at the terminal line run: 
#    % voms-proxy-init -voms cms
#    and enter your GRID passphrase. This is necessary to access the
#    miniAOD.
#
# 2. To run the script: % python plot_jet_properties.py
#
# Output pdf is in the same directory. 
#####################################################################
# Import necessary modules
import ROOT
import rootpy 
# import rootpy.plotting as plt
import pprint
from DataFormats.FWLite import Events, Handle
import pandas as pd 

import matplotlib
import numpy as numpy
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle   # for creating a custom legend
from matplotlib.ticker import AutoLocator, AutoMinorLocator, MultipleLocator, FormatStrFormatter # for ticks

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

#========= Load and format the data ==========# 

# Load the input file
data = []
events = Events(
    'root://xrootd-cms.infn.it//store/relval/CMSSW_9_2_7/RelValTTbarLepton_13/'
    'MINIAODSIM/PU25ns_92X_upgrade2017_realistic_v7-v1/00000/'
    '4E285340-6471-E711-AF0C-0025905B85F6.root'
)
# Create a handle (only needs to be done once, before the event loop)
handle = Handle('vector<pat::Jet>')
# List taggers you want to plot
taggers = [
	('CSV_AVR', 'pfCombinedSecondaryVertexV2BJetTags'),
	('deepCSV_BvsAll', 'pfDeepCSVDiscriminatorsJetTags:BvsAll'),
	('deepCSV_CvsB'  , 'pfDeepCSVDiscriminatorsJetTags:CvsB'),
	('deepCSV_CvsL'  , 'pfDeepCSVDiscriminatorsJetTags:CvsL'),
]

for event in events:
	event.getByLabel('slimmedJets', handle)
	jets = handle.product()
	for jet in jets:
		if jet.pt() < 20 or abs(jet.eta()) > 2.4: continue #basic selection
		entry = [
			abs(jet.hadronFlavour()), # 5 - b-jet, 4 - c-jet, 0 - light
			jet.pt(),
			jet.eta(),
			jet.phi(),
			jet.mass(),
		]
		entry.extend([jet.bDiscriminator(i) for _, i in taggers])
		data.append(entry)

# Make a DataFrame out of data
data = pd.DataFrame(data,
		    columns = ['flavour', 'pt', 'eta', 'phi', 'mass'] + 
		              [i for i, _ in taggers])
				   
# Create a 'isB' boolean tag to distinguish hadronFlavor == 5 and != 5 data
data['isB']  = (data.flavour == 5)
# Use the Pandas groupby method to create a dictionary with keys True and False.
# Now we can access all the b-jets with jets_dict[True], and all the non-b
# jets with jets_dict[False].
jets_dict = {k: v for k, v in data.groupby('isB')} 

#========= Plotting ===========================================#
### If adding more jet variables/subplots, something needs to be changed
### in the next 5 lines:
fig, axs = plt.subplots(4, 2,     # check no. of subplots 
			facecolor = 'w', edgecolor = 'k')

# State the variables we want to plot
plotFeatures = ['pt','eta','phi','mass'] + [i for i, _ in taggers] 
# x-axis range to plot
xRanges = [(0, 300), (-2.7, 2.7), (-3.5, 3.5), (0, 50),
	   (-0.1, 1.1), (0.1, 1.1), (0.1, 1.1), (0.1, 1.1)]

### Adjust font sizes, alpha, hatch style here
fsizeLarge = 12
fsizeTicks = 7  
alph = 0.8
htch = '.............'

# Return the axes as an array which we can access by 1 index
axs = axs.ravel() 

### Actual plotting starts here:
# Loop through the features and plot each one in a different subplot
for feature in plotFeatures:
	i = plotFeatures.index(feature)
	# In each subplot, plot non-bjet values in Blue and b-jet values in Red
	for isB, col in zip([False, True], ['blue', 'red']):
	        # Create subplot title
		axs[i].set_title(plotFeatures[i], fontsize = fsizeLarge)
		# Set x-axis range
		axs[i].set_xlim(xRanges[i])
		# Format axis ticks
		axs[i].xaxis.set_major_locator(AutoLocator())
		axs[i].yaxis.set_major_locator(AutoLocator())
		axs[i].xaxis.set_minor_locator(AutoMinorLocator())
		axs[i].yaxis.set_minor_locator(AutoMinorLocator())
		axs[i].tick_params(labelsize = fsizeTicks)
		# Create histogram from Pandas dataFrame and plot it
		jets_dict[isB][feature].hist(bins = 100,
					     ax = axs[i],
					     # Cosmetic options:
					     alpha = alph,
					     edgecolor = col,
					     fill = None,
					     grid = False,
					     hatch = htch,
					     histtype='step',
					     normed = True)

# Set common axes labels and title
fig.text(0.5, 0.01, 'Value', ha = 'center', va = 'center', 
	 fontsize = fsizeLarge)
fig.text(0.01, 0.5, 'Counts', ha = 'center', va = 'center', 
	 rotation = 'vertical', fontsize = fsizeLarge)
fig.text(0.5, 1, 'Jet properties in a simulated t-tbar event',
	 ha = 'center', va = 'center', fontsize = fsizeLarge)

### Create legend
handles = [Rectangle((0,0), 1, 1,
		     alpha = alph, edgecolor = c,
		     fill = None, hatch = htch) for c in ['blue', 'red']]
labels = ["Background (hadronFlavour != 5)", "Signal (hadronFlavour == 5)"]
# Switch current active axis to top right subplot
plt.axes(axs[1])
plt.legend(handles, labels,
	   # offset relative to top right subplot:
	   bbox_to_anchor = (1, 1.7),
	   fontsize = fsizeTicks,
	   loc = 'upper right')

# Adjust layout
plt.tight_layout()

# Save the figure and use bbox_inches to make room for legend
plt.savefig("jetPropertiesFromMiniAOD.pdf", bbox_inches = 'tight')
plt.clf()
		



