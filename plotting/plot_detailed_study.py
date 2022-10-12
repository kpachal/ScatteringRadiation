import ROOT
import math
import numpy as np
from art.morisot_2p0 import Morisot_2p0
import sys,os

# Initialize painter
myPainter = Morisot_2p0()
myPainter.colourpalette.setColourPalette("notSynthwave")
myPainter.labelType = 4
myPainter.savePDF = True
myPainter.luminosity = -1
myPainter.CME = -1

# Prediction for electron destinations (multiple scattering): goudsmit-saunderson

# The file to look at: just one
infilename = "../results_1micron_1e8events.root"

# For getting theory predictions: what foil thickness?
thickness = 1

# Where plots will go: default is subdir based on input file name
plotdir = "plots/{0}/".format(infilename.split("/")[-1].replace(".root",""))
#if not os.path.exists(plotdir) :
#  os.path.(plotdir)

# Only one particle present in input file? If so denote here
#useOnly = "neutron"
useOnly = ""

# Additional tag for output
tag = ""

# Neutron rebinning factor
rebin_neutrons = 5

def littleDrawStack(stack_list,data=None,signal_list=None,stack_labels=[],data_label="",xlabel="",ylabel="",plotname="",doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=2,extraLines=[],xl=None,xh=None,ylow=None,yhigh=None) :

    # Make a canvas
    c = myPainter.makeCanvas(plotname,logx,logy)

    goodcolours = myPainter.getFillColours(len(stack_list))

    # Make a stack for items in the stack_list
    stack = ROOT.THStack("stack","stacked histograms")
    for histogram in stack_list :
      index = stack_list.index(histogram)
      histogram.SetLineColor(goodcolours[index])
      histogram.SetLineWidth(2)
      histogram.SetFillColor(goodcolours[index])
      histogram.SetFillStyle(1001)
      histogram.SetTitle("")
      #histogram.GetYaxis().SetRangeUser(0.1,1e10)
      stack.Add(histogram,"hist")

    # Without data, stack sets plot ranges and formats
    stack.Draw("F")
    #stack.SetMaximum(stack.GetMaximum()*5.0 if logy else stack.GetMaximum()*1.4)
    stack.SetMinimum(0.5 if logy else 0)
    # A sensible stack should have the biggest contribution (i.e. longest tails) on top
    stack.GetXaxis().SetRangeUser(xl,xh)
    #stack.GetYaxis().SetRangeUser(0.1,stack.GetMaximum()*5.0)
    stack.GetXaxis().SetTitle(xlabel)
    stack.GetYaxis().SetTitle(ylabel)
    c.Update()

    # Finally, add some labels  
    c.SaveAs(plotname+".eps")
    c.SaveAs(plotname+".pdf")

def compute_rms_1D(hist) :
    sum2 = 0
    for bin in range(hist.GetNbinsX()+1) :
        sum2 += hist.GetBinContent(bin) * hist.GetBinCenter(bin)**2
    RMS = math.sqrt(sum2/hist.Integral())
    return RMS

def compute_rms_2D(hist) :
    sum2 = 0
    for binx in range(hist.GetNbinsX()+1) :
        for biny in range(hist.GetNbinsY()+1) :
          radius2 = hist.GetXaxis().GetBinCenter(binx)**2 + hist.GetYaxis().GetBinCenter(biny)**2
          sum2 += hist.GetBinContent(binx,biny) * radius2
    RMS = math.sqrt(sum2/hist.Integral())
    return RMS

savehists = {}
particle_dict = {"eminus" : {"full" : "electrons", "short" : "e-"},
                 "eplus" : {"full" : "positrons", "short" : "e+"},
                 "gamma" : {"full" : "photons", "short" : "#gamma"},
                 "neutron" : {"full" : "neutrons", "short" : "n"} }

# Pick up theory curves for e- multiple scattering
intheoryfile = ROOT.TFile.Open("../theory_curves.root","READ")
theory_curves = {}
for particle in ["eminus"] :
  thisdict = {}
  multiple_scattering_DCS = intheoryfile.Get("MDCS_{0}micron".format(thickness))
  thisdict["rad"] = multiple_scattering_DCS
  multiple_scattering_DCS_deg = intheoryfile.Get("MDCS_{0}micron_deg".format(thickness))
  thisdict["deg"] = multiple_scattering_DCS_deg
  theory_curves[particle] = thisdict
integrals_multiple_scattering = intheoryfile.Get("integrals_MDCS")
Widths_fromTDRFormula = intheoryfile.Get("Widths_fromTDRFormula")
Widths_fromTF1 = intheoryfile.Get("Widths_fromTF1")
intheoryfile.Close()

# Now get the main plots.

plot_dict = {}

infile = ROOT.TFile.Open(infilename)

for normalisation in [True, False] :

  for particle in ["eminus", "eplus", "gamma", "neutron"] :

    if (useOnly and particle not in useOnly) : 
      continue

    if normalisation :
      extratag = "_relative"
    else :
      extratag = ""

    particle_full = particle_dict[particle]["full"]
    particle_short = particle_dict[particle]["short"]

    # Normalised scattering angle
    polar_scattering = infile.Get("angle_polar_{0}".format(particle))
    polar_scattering.SetDirectory(0)
    if normalisation : polar_scattering.Scale(1.0/polar_scattering.GetMaximum())
    # For electrons: in what degree radius is 3.7 sigma contained?
    if particle == "eminus" :
      total = polar_scattering.Integral()
      sofar = 0
      for bin in range(1,polar_scattering.GetNbinsX()) :
        sofar += polar_scattering.GetBinContent(bin)
        if sofar/total > 0.986 :
          print("Upper edge is",polar_scattering.GetBinLowEdge(bin+1),"and percent is",sofar/total)
          break
    polar_scattering_norm = infile.Get("angle_polar_{0}_norm".format(particle))
    polar_scattering_norm.SetDirectory(0)
    if normalisation : polar_scattering_norm.Scale(1.0/polar_scattering_norm.GetMaximum())
    polar_scattering_norm.SetName(polar_scattering_norm.GetName()+"_{0}micron".format(thickness))

    # Normalised scattering angle, axis in deg
    polar_scattering_deg = infile.Get("angle_polar_deg_{0}".format(particle))
    polar_scattering_deg.SetDirectory(0)
    if normalisation : polar_scattering_deg.Scale(1.0/polar_scattering_deg.GetMaximum())
    # For electrons: in what degree radius is 3.7 sigma contained?
    if particle == "eminus" :
      total = polar_scattering_deg.Integral()
      sofar = 0
      for bin in range(1,polar_scattering_deg.GetNbinsX()) :
        sofar += polar_scattering_deg.GetBinContent(bin)
        if sofar/total > 0.986 :
          print("Upper edge is",polar_scattering_deg.GetBinLowEdge(bin+1),"and percent is",sofar/total)
          break    
    polar_scattering_deg_norm = infile.Get("angle_polar_deg_{0}_norm".format(particle))
    polar_scattering_deg_norm.SetDirectory(0)
    if normalisation : polar_scattering_deg_norm.Scale(1.0/polar_scattering_deg_norm.GetMaximum())
    polar_scattering_deg_norm.SetName(polar_scattering_deg_norm.GetName()+"_{0}micron".format(thickness))

    # Rebin neutrons: low stats,
    # but only after saving raw version so plots stack better
    # if ("neutron" in particle) :
    #     polar_scattering.Rebin(rebin_neutrons)
    #     polar_scattering_norm.Rebin(rebin_neutrons)
    #     polar_scattering_deg.Rebin(rebin_neutrons)
    #     polar_scattering_deg_norm.Rebin(rebin_neutrons)    

    # Save
    plot_dict["polar_scattering_{0}".format(particle)] = polar_scattering
    plot_dict["polar_scattering_norm_{0}".format(particle)] = polar_scattering_norm    
    plot_dict["polar_scattering_deg_{0}".format(particle)] = polar_scattering_deg
    plot_dict["polar_scattering_deg_norm_{0}".format(particle)] = polar_scattering_deg_norm

    # Make individual plots by thickness
    # Radians, unnormalised
    yhigh = None
    if "neutron" in particle :
        yhigh = 10

    myPainter.drawOverlaidHistos([polar_scattering],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{1}umber of {0} [1/rad]".format(particle_full,("Relative n" if normalisation else "N")),plotname=plotdir+"scattering_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=3.14,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{1}umber of {0} [1/rad]".format(particle_full,("Relative n" if normalisation else "N")),plotname=plotdir+"scattering_{0}_{1}microns{2}_zoom".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1.0,ylow=None,yhigh=yhigh,useTrueEdges=True)    
    # Radians, normalised
    myPainter.drawOverlaidHistos([polar_scattering_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{1}umber of {0} [1/str]".format(particle_full,("Relative n" if normalisation else "N")),plotname=plotdir+"scattering_norm_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=3.14,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{1}umber of {0} [1/str]".format(particle_full,("Relative n" if normalisation else "N")),plotname=plotdir+"scattering_norm_{0}_{1}microns{2}_zoom".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1.0,ylow=None,yhigh=yhigh,useTrueEdges=True)
    # Degrees, unnormalised
    myPainter.drawOverlaidHistos([polar_scattering_deg],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} number of {1}".format(("Relative" if normalisation else "Total"), particle_full),plotname=plotdir+"scattering_deg_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=180,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_deg],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} number of {1}".format(("Relative" if normalisation else "Total"), particle_full),plotname=plotdir+"scattering_deg_{0}_{1}microns{2}_zoom".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=40,ylow=None,yhigh=yhigh,useTrueEdges=True)   
    # Degrees, normalised
    myPainter.drawOverlaidHistos([polar_scattering_deg_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{1}umber of {0} [1/deg^{2}]".format(particle_full,("Relative n" if normalisation else "N"),"{2}"),plotname=plotdir+"scattering_deg_norm_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=180,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_deg_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{1}umber of {0} [1/deg^{2}]".format(particle_full,("Relative n" if normalisation else "N"),"{2}"),plotname=plotdir+"scattering_deg_norm_{0}_{1}microns{2}_zoom".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=40,ylow=None,yhigh=yhigh,useTrueEdges=True) 

    # Now normalize to match theory curves and plot.
    # Only this if e-
    # Currently only good for non-normalised. examine for otherwise.
    if "eminus" in particle and not normalisation :
        # Should already be fitted to the distribution
        rad_DCS = theory_curves[particle]["rad"]
        myPainter.drawHistsWithTF1s([polar_scattering_norm],[rad_DCS],as_data=False,match_colours=True,hist_labels=["Geant4 e-"],func_labels=["Moliere DCS"],xlabel="Polar scattering angle [rad.]",ylabel="Arbitrary units",plotname=plotdir+"scattering_w_theory_{0}microns".format(thickness)+tag,doRatios=False,ratioName="",doErrors=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1.0,ylow=None,yhigh=1e3)
        deg_DCS = theory_curves[particle]["deg"]
        myPainter.drawHistsWithTF1s([polar_scattering_deg_norm],[deg_DCS],as_data=False,match_colours=True,hist_labels=["Geant4 e-"],func_labels=["Moliere DCS"],xlabel="Polar scattering angle [deg.]",ylabel="Arbitrary units",plotname=plotdir+"scattering_w_theory_deg_{0}microns".format(thickness)+tag,doRatios=False,ratioName="",doErrors=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=90.,ylow=None,yhigh=1e3)            

    # Save positions and momenta
    for histname in ["momentum_x_{0}".format(particle), "momentum_y_{0}".format(particle), "momentum_z_{0}".format(particle), "position_x_{0}".format(particle), "position_y_{0}".format(particle), "position_z_{0}".format(particle)] :
        thishist = infile.Get(histname)
        thishist.SetDirectory(0)
        if normalisation : thishist.Scale(1.0/thishist.GetMaximum())
        plot_dict[histname] = thishist

    # p
    momentum = infile.Get("momentum_{0}".format(particle))
    momentum.SetDirectory(0)
    if normalisation : momentum.Scale(1.0/momentum.GetMaximum())
    #if "neutron" in particle : momentum.Rebin(rebin_neutrons)
    plot_dict["momentum_{0}".format(particle)] = momentum
    myPainter.drawOverlaidHistos([momentum],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Momentum [MeV]",ylabel="{1}umber of {0}".format(particle_full,("Relative n" if normalisation else "N")),plotname=plotdir+"momentum_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=None,ylow=None,yhigh=None,useTrueEdges=True)

    # Energy (not representative for neutrons)
    energy = infile.Get("energy_{0}".format(particle))
    energy.SetDirectory(0)
    if normalisation : momentum.Scale(1.0/momentum.GetMaximum())
    plot_dict["energy_{0}".format(particle)] = energy
    myPainter.drawOverlaidHistos([energy],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Energy [MeV]",ylabel="{1}umber of {0}".format(particle_full,("Relative n" if normalisation else "N")),plotname=plotdir+"energy_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=None,ylow=None,yhigh=None,useTrueEdges=True)

    # Collect and save 2d histograms
    momentum_2d = infile.Get("momentum_xy_{0}".format(particle))
    momentum_2d.SetDirectory(0)
    if normalisation : momentum_2d.Scale(1.0/momentum_2d.GetMaximum())
    plot_dict["momentum_xy_{0}".format(particle)] = momentum_2d
    position_2d = infile.Get("position_xy_{0}".format(particle))
    position_2d.SetDirectory(0)
    if normalisation : position_2d.Scale(1.0/position_2d.GetMaximum())
    plot_dict["position_xy_{0}".format(particle)] = position_2d

    # Plot for visual fun
    myPainter.draw2DHist(momentum_2d,plotdir+"momentum_xy_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,"px [MeV]",-30,30,"py [MeV]",-30,30,particle_full,logz=True)
    myPainter.draw2DHist(position_2d,plotdir+"position_xy_{0}_{1}microns{2}".format(particle,thickness,extratag)+tag,"x [mm]",-1000,1000,"y [mm]",-1000,1000,particle_full,logz=True)

  savehists[particle] = plot_dict

# Make stacks of output particles within thickness. Only if doing all particles
if not useOnly :
  to_stack = [plot_dict["polar_scattering_norm_neutron"],plot_dict["polar_scattering_norm_eplus"],plot_dict["polar_scattering_norm_gamma"],plot_dict["polar_scattering_norm_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [rad]",ylabel="Norm. particles per rad^{2}",plotname=plotdir+"stacked_particles_norm_{0}microns".format(thickness)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1,ylow=None,yhigh=None)

  to_stack = [plot_dict["polar_scattering_neutron"],plot_dict["polar_scattering_eplus"],plot_dict["polar_scattering_gamma"],plot_dict["polar_scattering_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [rad]",ylabel="Total particles",plotname=plotdir+"stacked_particles_{0}microns".format(thickness)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1,ylow=None,yhigh=None)

  to_stack = [plot_dict["polar_scattering_deg_norm_neutron"],plot_dict["polar_scattering_deg_norm_eplus"],plot_dict["polar_scattering_deg_norm_gamma"],plot_dict["polar_scattering_deg_norm_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [deg]",ylabel="Norm. particles per deg^{2}",plotname=plotdir+"stacked_particles_deg_norm_{0}microns".format(thickness)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=90,ylow=None,yhigh=None)

  to_stack = [plot_dict["polar_scattering_deg_neutron"],plot_dict["polar_scattering_deg_eplus"],plot_dict["polar_scattering_deg_gamma"],plot_dict["polar_scattering_deg_eminus"]]
  littleDrawStack(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [deg]",ylabel="Total particles",plotname=plotdir+"stacked_particles_deg_{0}microns".format(thickness)+tag,doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xl=0,xh=90,ylow=None,yhigh=None)

# Close input file
infile.Close()

# Save histograms that I created here, so I can look at them further if desired
outfile = ROOT.TFile.Open("plots_detailedstudy.root","RECREATE")
outfile.cd()
for thickness in savehists.keys() :
  for name,hist in savehists[thickness].items() :
    hist.Write()
outfile.Close()