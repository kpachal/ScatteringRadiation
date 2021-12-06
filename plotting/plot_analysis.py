import ROOT
import math
import numpy as np
from art.morisot_2p0 import Morisot_2p0

# Initialize painter
myPainter = Morisot_2p0()
myPainter.colourpalette.setColourPalette("notSynthwave")
myPainter.labelType = 4
myPainter.savePDF = True
# 0 Just ATLAS    
# 1 "Preliminary"
# 2 "Internal"
# 3 "Simulation Preliminary"
# 4 "Simulation Internal"
# 5 "Simulation"
# 6 "Work in Progress"

# Prediction for electron destinations (multiple scattering): goudsmit-saunderson

infilename = "../results_{0}micron_1e8events.root"

def compute_rms_1D(hist) :
    sum2 = 0
    for bin in range(hist.GetNbinsX()+1) :
        sum2 += hist.GetBinContent(bin) * hist.GetBinCenter(bin)**2
    RMS = math.sqrt(sum2/hist.Integral())
    return RMS

def normalise_solid_angle(inhist, radians=True) :

    newhist = inhist.Clone(inhist.GetName()+"_norm")
    newhist.SetDirectory(0)
    newhist.Reset()

    # Solid angle subtended = 2 pi [-cos theta] | theta 1 to theta 2
    for bin in range(newhist.GetNbinsX()+1) :
        fullval = inhist.GetBinContent(bin)
        theta1 = inhist.GetBinLowEdge(bin)
        theta2 = inhist.GetBinLowEdge(bin+1)
        if radians :
            normval = 2 * math.pi * (- math.cos(theta2) + math.cos(theta1))
        else :
            normval = 2 * math.pi * (- math.cos(np.radians(theta2)) + math.cos(np.radians(theta1)))
        newhist.SetBinContent(bin,fullval/normval)
        newhist.SetBinError(bin,inhist.GetBinError(bin)/normval)

    return newhist

savehists = {}
particle_dict = {"eminus" : {"full" : "Electrons", "short" : "e-"},
                 "eplus" : {"full" : "Positrons", "short" : "e+"},
                 "gamma" : {"full" : "Photons", "short" : "#gamma"},
                 "neutron" : {"full" : "Neutrons", "short" : "n"} }
for thickness in 1, 5, 10 :

  thisThickness = {}

  thisname = infilename.format(thickness)
  print("Opening file",thisname)
  infile = ROOT.TFile.Open(thisname)

  for particle in ["eminus", "eplus", "gamma", "neutron"] :

    particle_full = particle_dict[particle]["full"]
    particle_short = particle_dict[particle]["short"]

    # Normalised scattering angle
    polar_scattering = infile.Get("angle_polar_{0}".format(particle))
    polar_scattering.SetDirectory(0)
    polar_scattering_norm = normalise_solid_angle(polar_scattering,radians=True)
    polar_scattering_norm.SetName(polar_scattering_norm.GetName()+"_{0}micron".format(thickness))
    thisThickness["polar_scattering_{0}".format(particle)] = polar_scattering
    thisThickness["polar_scattering_norm_{0}".format(particle)] = polar_scattering_norm

    # Normalised scattering angle, axis in deg
    polar_scattering_deg = infile.Get("angle_polar_deg_{0}".format(particle))
    polar_scattering_deg.SetDirectory(0)
    polar_scattering_deg_norm = normalise_solid_angle(polar_scattering_deg,radians=False)
    polar_scattering_deg_norm.SetName(polar_scattering_deg_norm.GetName()+"_{0}micron".format(thickness))
    thisThickness["polar_scattering_deg_{0}".format(particle)] = polar_scattering_deg
    thisThickness["polar_scattering_deg_norm_{0}".format(particle)] = polar_scattering_deg_norm

    # Make individual plots by thickness
    # Radians, unnormalised
    yHigh = None
    if "neutron" in particle :
        yHigh = 10

    myPainter.drawOverlaidHistos([polar_scattering],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} total".format(particle_full),plotname="scattering_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=3.14,yLow=None,yHigh=yHigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} total".format(particle_full),plotname="scattering_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=1.0,yLow=None,yHigh=yHigh,useTrueEdges=True)    
    # Radians, normalised
    myPainter.drawOverlaidHistos([polar_scattering_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_norm_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=3.14,yLow=None,yHigh=yHigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_norm_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=1.0,yLow=None,yHigh=yHigh,useTrueEdges=True)
    # Degrees, unnormalised
    myPainter.drawOverlaidHistos([polar_scattering_deg],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} total".format(particle_full),plotname="scattering_deg_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=180,yLow=None,yHigh=yHigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_deg],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} total".format(particle_full),plotname="scattering_deg_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=40,yLow=None,yHigh=yHigh,useTrueEdges=True)   
    # Degrees, normalised
    myPainter.drawOverlaidHistos([polar_scattering_deg_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_deg_norm_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=180,yLow=None,yHigh=yHigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_deg_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_deg_norm_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=40,yLow=None,yHigh=yHigh,useTrueEdges=True) 

    # p (energy not representative for neutrons)
    momentum = infile.Get("momentum_{0}".format(particle))
    momentum.SetDirectory(0)
    thisThickness["momentum_{0}".format(particle)] = momentum
    myPainter.drawOverlaidHistos([momentum],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Momentum [MeV]",ylabel="Number of {0}".format(particle_full),plotname="momentum_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=None,yLow=None,yHigh=None,useTrueEdges=True)

  savehists[thickness] = thisThickness

  # Make stacks of output particles within thickness
  to_stack = [thisThickness["polar_scattering_norm_neutron"],thisThickness["polar_scattering_norm_eplus"],thisThickness["polar_scattering_norm_gamma"],thisThickness["polar_scattering_norm_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [rad]",ylabel="Norm. particles per rad^{2}",plotname="stacked_particles_norm_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=1,yLow=None,yHigh=None)

  to_stack = [thisThickness["polar_scattering_neutron"],thisThickness["polar_scattering_eplus"],thisThickness["polar_scattering_gamma"],thisThickness["polar_scattering_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [rad]",ylabel="Total particles",plotname="stacked_particles_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=1,yLow=None,yHigh=None)

  to_stack = [thisThickness["polar_scattering_deg_norm_neutron"],thisThickness["polar_scattering_deg_norm_eplus"],thisThickness["polar_scattering_deg_norm_gamma"],thisThickness["polar_scattering_deg_norm_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [deg]",ylabel="Norm. particles per deg^{2}",plotname="stacked_particles_deg_norm_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=90,yLow=None,yHigh=None)

  to_stack = [thisThickness["polar_scattering_deg_neutron"],thisThickness["polar_scattering_deg_eplus"],thisThickness["polar_scattering_deg_gamma"],thisThickness["polar_scattering_deg_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [deg]",ylabel="Total particles",plotname="stacked_particles_deg_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=90,yLow=None,yHigh=None)

# Make plots comparing thicknesses
histos_list = [savehists[1]["polar_scattering_deg_eminus"],savehists[5]["polar_scattering_deg_eminus"],savehists[10]["polar_scattering_deg_eminus"]]
myPainter.drawOverlaidHistos(histos_list,data=None,signal_list=None,histos_labels=["1 #mum foil","5 #mum foil","10 #mum foil"],data_label="",xlabel="Polar angle [deg]",ylabel="Total electrons",plotname="electron_scattering_deg_compareThickness".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=45,yLow=None,yHigh=None)
histos_list = [savehists[1]["polar_scattering_deg_norm_eminus"],savehists[5]["polar_scattering_deg_norm_eminus"],savehists[10]["polar_scattering_deg_norm_eminus"]]
myPainter.drawOverlaidHistos(histos_list,data=None,signal_list=None,histos_labels=["1 #mum foil","5 #mum foil","10 #mum foil"],data_label="",xlabel="Polar angle [deg]",ylabel="Norm. electrons per deg^{2}",plotname="electron_scattering_deg_norm_compareThickness".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,luminosity=-1,CME=-1,nLegendColumns=1,extraLines=[],xLow=0,xHigh=45,yLow=None,yHigh=None)

# Save histograms that I created here, so I can look at them further if desired
outfile = ROOT.TFile.Open("plots.root","RECREATE")
outfile.cd()
for thickness in savehists.keys() :
  for name,hist in savehists[thickness].items() :
    hist.Write()
outfile.Close()

# Beam spread results, electrons only - compare methods
# Not using ROOT's RMS since it compares to the mean along the x axis and I want to compare to zero
for thickness in 1, 5, 10 :
  # RMS of normalised distribution
  RMS_norm_dist = compute_rms_1D(savehists[thickness]["polar_scattering_norm_eminus"])
  print("RMS from solid-angle-normalised distribution", RMS_norm_dist,"=",math.degrees(RMS_norm_dist),"degrees")
  # RMS by hand using bin center - non-normalised distribution
  RMS_notnormalised = compute_rms_1D(savehists[thickness]["polar_scattering_eminus"])
  print("RMS from unnormalised distribution", RMS_notnormalised,"=",math.degrees(RMS_notnormalised),"degrees")
  # RMS of 2D distribution - position - and calculate from there


  #RMS_norm_dist = compute_rms_1D(savehists[thickness]["polar_scattering_deg_norm_eminus"])
#   # RMS by hand using bin center - non-normalised distribution
#   RMS_notnormalised = compute_rms_1D(savehists[thickness]["polar_scattering_eminus"])
#   print("RMS from full distribution is",math.degrees(RMS_notnormalised))
#   print("RMS from solid-angle-normalised distribution is",math.degrees(RMS_norm_dist))