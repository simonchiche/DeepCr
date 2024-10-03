# This script reads the 2D NumberOfParticlesDepth histos from a g4 simulation and plots
# the number of e+/e- in function of depth

import sys
import ROOT

if(len(sys.argv) != 3):
    sys.exit("This script expects 2 arguments:\n - The root file containing the histograms.\n - The output directory to store the plot in.")
else:
    root_file = sys.argv[1]
    output_dir = sys.argv[2]

# Reading the root file
input_file = ROOT.TFile(root_file, "OLD")

# The histos_combined.root file has 2D histograms indicating the number of particles in the in-ice shower
# with depth along the shower axis (in m) on the x-axis 
# and kinetic energy (in MeV) on the y-axis.
# There are six such histograms:
# - NumberOfParticlesDepth_gammas, which indicates the number of gammas in the in-ice shower
# - NumberOfParticlesDepth_electrons, which indicates the number of electrons in the in-ice shower
# - NumberOfParticlesDepth_positrons, which indicates the number of positrons in the in-ice shower
# - NumberOfParticlesDepth_muons, which indicates the number of muons in the in-ice shower
# - NumberOfParticlesDepth_antimuons, which indicates the number of antimuons in the in-ice shower
# - NumberOfParticlesDepth_hadrons, which indicates the number of hadrons in the in-ice shower
# Summing these histograms together indicates the total number of particles in the in-ice shower.
# Here we're interested in the total electrons + positrons:
nop_tot_2D = input_file.Get("NumberOfParticlesDepth_positrons")
nop_tot_2D.Add(input_file.Get("NumberOfParticlesDepth_electrons"))

# Use these to add more types of particles:
#nop_tot_2D.Add(input_file.Get("NumberOfParticlesDepth_gammas"))
#nop_tot_2D.Add(input_file.Get("NumberOfParticlesDepth_antimuons"))
#nop_tot_2D.Add(input_file.Get("NumberOfParticlesDepth_muons"))
#nop_tot_2D.Add(input_file.Get("NumberOfParticlesDepth_hadrons"))

# Giving the new histogram a proper title, just for bookkeeping
nop_tot_2D.SetTitle("Number of e+/e- in function of depth")

# Doing the integration over the y-axis (energy), so we end up with a 1D histogram simply indicating
# number of particles in function of depth (in m) along the shower axis.
nop_tot_1D = nop_tot_2D.ProjectionX()

# Getting the values corresponding to the maximum number of particles (depth, number of particles and 
# bin number).
# Note the unusual bin numbering scheme of ROOT:
# - bins run from "1" to "number of bins", instead the usual "0" to "number of bins - 1"
# - bins "0" and "number of bins + 1" hold the overflow, i.e. all the particles outside the xrange
depthval_max = -999.
nopval_max = -999.
binnr_max = -999
for i in range(1, nop_tot_1D.GetXaxis().GetNbins()+1):
    if(nop_tot_1D.GetBinContent(i) > nopval_max):
        nopval_max = nop_tot_1D.GetBinContent(i)
        depthval_max = nop_tot_1D.GetBinCenter(i)
        binnr_max = i

print("Xmax (slant) is %s m" % depthval_max)

# Preparing to plot the histogram
ROOT.gStyle.SetOptStat(0) # Removes the statistics box in the plot
ROOT.gROOT.SetBatch(ROOT.kTRUE) # Disables display pop-up of canvas
ROOT.gStyle.SetTitleXSize(0.03)
ROOT.gStyle.SetTitleYSize(0.03)

# Copying some settings from the original histograms in the histos_combined.root file
nop_tot_1D.GetXaxis().SetTitleFont(nop_tot_1D.GetYaxis().GetTitleFont())
nop_tot_1D.GetXaxis().SetLabelFont(nop_tot_1D.GetYaxis().GetLabelFont())
nop_tot_1D.GetXaxis().SetLabelSize(nop_tot_1D.GetYaxis().GetLabelSize())

# Line width of the plot
nop_tot_1D.SetLineWidth(2)

# The canvas which holds the plot
cv = ROOT.TCanvas("cv", "cv", 1500, 900)
cv = ROOT.TCanvas("cv", "cv", 1500, 900)
cv.SetLeftMargin(0.34)
cv.SetRightMargin(0.34)
cv.SetTopMargin(0.24)
cv.SetBottomMargin(0.24)

# Making and saving the plot
nop_tot_1D.Draw("HIST")
cv.SaveAs("%s/nop_depth.png" % output_dir)

# Closing the input file
input_file.Close()
