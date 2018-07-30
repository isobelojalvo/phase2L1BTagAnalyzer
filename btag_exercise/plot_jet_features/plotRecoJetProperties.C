//// This script plots reconstructed jet properties from a TTree
//// file created by the phase2L1BTagTrigger analyzer. 
////
//// Usage: 
//// % root -l
//// % .x plotRecoJetProperties.C
////  
//// Output .PDF file is in the same directory.
////
//// Based on mlpHiggs.C in official ROOT tutorials.
//// Takes several minutes to run on one .ROOT file.

void plotRecoJetProperties() {
	////////////////////////////
	/// Load data and initialize
	////////////////////////////

	const char *fname = "/afs/cern.ch/work/s/skkwan/public/triggerDevel/CMSSW_10_1_5/src/L1Trigger/phase2L1BTagAnalyzer/test/analyzer.root";
	// Load file with checks on whether it exists
	TFile *input = 0;
	if (!gSystem->AccessPathName(fname)) {
		input = TFile::Open(fname);
	} else {
		printf("accessing %s file",fname);
		input = TFile::Open(Form("%s", fname));
	}
	if (!input) return;

	// Read in data
	TTree *raw = (TTree *) input->Get("L1BTagAnalyzer/efficiencyTree");

	// New TTrees for sig and background 
	TTree *sig = new TTree("sig", "hadronFlavor == 5");
	TTree *bkg = new TTree("bkg", "hadronFlavor != 5");

	// Declare variables to read from TTree
	Double_t recoPt, recoEta, recoPhi, recoTk1IP, recoTk2IP, recoTk3IP, recoTk4IP, l1Pt, l1Eta, l1Phi;
	Int_t hadronFlavor;
   
	// Set branch addresses
	raw->SetBranchAddress("recoPt",  &recoPt);
	raw->SetBranchAddress("recoEta", &recoEta);
	raw->SetBranchAddress("recoPhi", &recoPhi);
	raw->SetBranchAddress("hadronFlavor", &hadronFlavor);
	raw->SetBranchAddress("recoTk1IP", &recoTk1IP);
	raw->SetBranchAddress("recoTk2IP", &recoTk2IP);
	raw->SetBranchAddress("recoTk3IP", &recoTk3IP);
	raw->SetBranchAddress("recoTk4IP", &recoTk4IP);

	//////////////////////////////////////////////////////////////////////
	// Set up histograms
	// note the x-axis range is manually set: best to check binning and
	// ranges in TBrowser / interactively in ROOT as a first pass.
	//////////////////////////////////////////////////////////////////////

	// Global options (see TStyle documentation for all options)
	gStyle->SetOptStat(0); 

	// Signal histograms
	gStyle->SetHistLineColor(kRed); 
	gStyle->SetHistFillColor(kRed);
	gStyle->SetHistFillStyle(3001);	 
	TH1D *sig_recoPt = new TH1D("sig_recoPt", "sig_recoPt", 100, 0, 175);
	TH1D *sig_recoEta = new TH1D("sig_recoEta", "sig_recoEta", 100, -2.7, 2.7);
	TH1D *sig_recoPhi = new TH1D("sig_recoPhi", "sig_recoPhi", 100, -3.1, 3.1);
	TH1D *sig_recoTk1IP = new TH1D("sig_recoTk1IP", "sig_recoTk1IP", 70, -0.1, 0.1);
	TH1D *sig_recoTk2IP = new TH1D("sig_recoTk2IP", "sig_recoTk2IP", 70, -0.5, 0.5);
	TH1D *sig_recoTk3IP = new TH1D("sig_recoTk3IP", "sig_recoTk3IP", 700, -0.5, 0.5);
	TH1D *sig_recoTk4IP = new TH1D("sig_recoTk4IP", "sig_recoTk4IP", 70, -0.5, 0.5);
		
	// Likewise, background histograms
	gStyle->SetHistLineColor(kBlue); 
	gStyle->SetHistFillColor(kBlue);
	gStyle->SetHistFillStyle(3001);
	TH1D *bkg_recoPt = new TH1D("bkg_recoPt", "bkg recoPt", 100, 0, 175);
	TH1D *bkg_recoEta = new TH1D("bkg_recoEta", "bkg_recoEta", 100, -2.7, 2.7);
	TH1D *bkg_recoPhi = new TH1D("bkg_recoPhi", "bkg_recoPhi", 100, -3.1, 3.1);
	TH1D *bkg_recoTk1IP = new TH1D("bkg_recoTk1IP", "sig_recoTk1IP", 70, -0.1, 0.1);
	TH1D *bkg_recoTk2IP = new TH1D("bkg_recoTk2IP", "bkg_recoTk2IP", 70, -0.5, 0.5);
	TH1D *bkg_recoTk3IP = new TH1D("bkg_recoTk3IP", "bkg_recoTk3IP", 70, -0.5, 0.5);
	TH1D *bkg_recoTk4IP = new TH1D("bkg_recoTk4IP", "bkg_recoTk4IP", 70, -0.5, 0.5);

	// Loop through jets and use hadronFlavor to distinguish signal and background
	Int_t i;
	for (i = 0; i < raw->GetEntries(); i++) {
		raw->GetEntry(i);
		if (hadronFlavor == 5) { // Signal
			sig_recoPt->Fill(recoPt);
			sig_recoEta->Fill(recoEta);
			sig_recoPhi->Fill(recoPhi);
			sig_recoTk1IP->Fill(recoTk1IP);
			sig_recoTk2IP->Fill(recoTk2IP);
			sig_recoTk3IP->Fill(recoTk3IP);
			sig_recoTk4IP->Fill(recoTk4IP);
		}
		else {	// background
			bkg_recoPt ->Fill(recoPt);
			bkg_recoEta->Fill(recoEta);
			bkg_recoPhi->Fill(recoPhi);
			bkg_recoTk1IP->Fill(recoTk1IP);
			bkg_recoTk2IP->Fill(recoTk2IP);
			bkg_recoTk3IP->Fill(recoTk3IP);
			bkg_recoTk4IP->Fill(recoTk4IP);
		}
	}

	//////////////////////////
	/// Plotting
	//////////////////////////
	// Create and set up canvas for plotting sig vs. bkg variables
	TCanvas* sb_canvas = new TCanvas("sb_canvas", "Reconstructed jet variables");
	sb_canvas->Divide(2, 4);
	

	sb_canvas->cd(1); sb_canvas->SetTitle("recoPt");
	sig_recoPt->DrawNormalized();  bkg_recoPt->DrawNormalized("same");

	sb_canvas->cd(2); sb_canvas->SetTitle("reco Eta");
	sig_recoEta->DrawNormalized(); bkg_recoEta->DrawNormalized("same"); 

	sb_canvas->cd(3); sb_canvas->SetTitle("reco Phi");
	sig_recoPhi->DrawNormalized(); bkg_recoPhi->DrawNormalized("same"); 

	sb_canvas->cd(4); sb_canvas->SetTitle("recoTk1IP");
	sig_recoTk1IP->DrawNormalized(); bkg_recoTk1IP->DrawNormalized("same"); 

	sb_canvas->cd(5); sb_canvas->SetTitle("recoTk2IP");
	sig_recoTk2IP->DrawNormalized(); bkg_recoTk2IP->DrawNormalized("same"); 

	sb_canvas->cd(6); sb_canvas->SetTitle("recoTk3IP");
	sig_recoTk3IP->DrawNormalized(); bkg_recoTk3IP->DrawNormalized("same"); 

	sb_canvas->cd(7); sb_canvas->SetTitle("recoTk4IP");
	sig_recoTk4IP->DrawNormalized(); bkg_recoTk4IP->DrawNormalized("same"); 

	// Create legend
	sb_canvas->cd(8);
	TLegend *legend = new TLegend(0.2, 0.2, .6, .6);
	legend->AddEntry(sig_recoPt, "Signal (hadronFlavor == 5)");
	legend->AddEntry(bkg_recoPt, "Background (hadronFlavor != 5)");
	legend->Draw();

	sb_canvas->SaveAs("output_plotRecoJetProperties.pdf");

	input->Close();
}
