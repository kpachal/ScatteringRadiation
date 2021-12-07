import ROOT
import math
import numpy as np
from art.morisot_2p0 import Morisot_2p0

# Initialize painter
myPainter = Morisot_2p0()
myPainter.colourpalette.setColourPalette("notSynthwave")
myPainter.labelType = 4
myPainter.savePDF = True
myPainter.luminosity = -1
myPainter.CME = -1
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

def compute_rms_2D(hist) :
    sum2 = 0
    for binx in range(hist.GetNbinsX()+1) :
        for biny in range(hist.GetNbinsY()+1) :
          radius2 = hist.GetXaxis().GetBinCenter(binx)**2 + hist.GetYaxis().GetBinCenter(biny)**2
          sum2 += hist.GetBinContent(binx,biny) * radius2
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

def norm_int_by_binwidth(hist) :

    sum = 0
    for bin in range(hist.GetNbinsX()+1) :
        sum = sum + hist.GetBinContent(bin)/hist.GetBinWidth(bin)
    return sum

savehists = {}
particle_dict = {"eminus" : {"full" : "Electrons", "short" : "e-"},
                 "eplus" : {"full" : "Positrons", "short" : "e+"},
                 "gamma" : {"full" : "Photons", "short" : "#gamma"},
                 "neutron" : {"full" : "Neutrons", "short" : "n"} }

# Pick up theory curves for e- multiple scattering
intheoryfile = ROOT.TFile.Open("../theory_curves.root","READ")
theory_curves = {}
for thickness in 1, 5, 10 :
  multiple_scattering_DCS = intheoryfile.Get("MDCS_{0}micron".format(thickness))
  theory_curves[thickness] = multiple_scattering_DCS
integrals_multiple_scattering = intheoryfile.Get("integrals_MDCS")
intheoryfile.Close()

# Now get the main plots.
for thickness in 1, 5, 10 :

  thisThickness = {}
  iThickness = [1, 5, 10].index(thickness)

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
    yhigh = None
    if "neutron" in particle :
        yhigh = 10

    myPainter.drawOverlaidHistos([polar_scattering],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} total".format(particle_full),plotname="scattering_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=3.14,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} total".format(particle_full),plotname="scattering_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1.0,ylow=None,yhigh=yhigh,useTrueEdges=True)    
    # Radians, normalised
    myPainter.drawOverlaidHistos([polar_scattering_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_norm_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=3.14,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [rad]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_norm_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1.0,ylow=None,yhigh=yhigh,useTrueEdges=True)
    # Degrees, unnormalised
    myPainter.drawOverlaidHistos([polar_scattering_deg],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} total".format(particle_full),plotname="scattering_deg_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=180,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_deg],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} total".format(particle_full),plotname="scattering_deg_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=40,ylow=None,yhigh=yhigh,useTrueEdges=True)   
    # Degrees, normalised
    myPainter.drawOverlaidHistos([polar_scattering_deg_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_deg_norm_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=180,ylow=None,yhigh=yhigh,useTrueEdges=True)
    myPainter.drawOverlaidHistos([polar_scattering_deg_norm],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Polar angle [deg]",ylabel="{0} per unit solid angle".format(particle_full),plotname="scattering_deg_norm_{0}_{1}microns_zoom".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=40,ylow=None,yhigh=yhigh,useTrueEdges=True) 

    # Now normalize to match theory curves and plot.
    # Only this if e-
    if "eminus" in particle :
        this_DCS = theory_curves[thickness]
        # Believe it should match solid-angle-normalized.
        # Match to bin edges for safety.
        startEdge = polar_scattering_norm.GetBinLowEdge(polar_scattering_norm.FindBin(0.01)+1)
        stopEdge = polar_scattering_norm.GetBinLowEdge(polar_scattering_norm.FindBin(1.0))
        integral_theory = integrals_multiple_scattering[iThickness]
        integral_current = polar_scattering_norm.Integral(polar_scattering_norm.FindBin(0.01),polar_scattering_norm.FindBin(1.0))
        print("Integral from bin",polar_scattering_norm.FindBin(0.01),"to bin",polar_scattering_norm.FindBin(1.0),"is",integral_current)
        print("whereas total integral is",polar_scattering_norm.Integral(1,polar_scattering_norm.GetNbinsX()))
        print("Integral norm to bin width is",norm_int_by_binwidth(polar_scattering_norm))
        print("Want to scale to",integral_theory)
        polar_scattering_norm.Scale(integral_theory/integral_current)
        myPainter.drawHistsWithTF1s([polar_scattering_norm],[this_DCS],as_data=False,match_colours=True,hist_labels=["Geant4 e-"],func_labels=["Moliere DCS"],xlabel="Polar scattering angle [rad.]",ylabel="Arbitrary units",plotname="scattering_w_theory_{0}microns".format(thickness),doRatios=False,ratioName="",doErrors=False,logx=False,logy=True,nLegendColumns=2,extraLines=[],xlow=0,xhigh=1.0,ylow=None,yhigh=1e3)    

    # Save positions and momenta
    momx = infile.Get("momentum_x_{0}".format(particle))
    momx.SetDirectory(0)
    thisThickness["momentum_x_{0}".format(particle)] = momx
    momy = infile.Get("momentum_y_{0}".format(particle))
    momy.SetDirectory(0)
    thisThickness["momentum_y_{0}".format(particle)] = momy
    momz = infile.Get("momentum_z_{0}".format(particle))
    momz.SetDirectory(0)
    thisThickness["momentum_z_{0}".format(particle)] = momz
    posx = infile.Get("position_x_{0}".format(particle))
    posx.SetDirectory(0)
    thisThickness["position_x_{0}".format(particle)] = posx
    posy = infile.Get("position_y_{0}".format(particle))
    posy.SetDirectory(0)
    thisThickness["position_y_{0}".format(particle)] = posy
    posz = infile.Get("position_z_{0}".format(particle))
    posz.SetDirectory(0)
    thisThickness["position_z_{0}".format(particle)] = posz

    # p (energy not representative for neutrons)
    momentum = infile.Get("momentum_{0}".format(particle))
    momentum.SetDirectory(0)
    thisThickness["momentum_{0}".format(particle)] = momentum
    myPainter.drawOverlaidHistos([momentum],data=None,signal_list=None,histos_labels=["Geant4 {0}".format(particle_short)],data_label="",xlabel="Momentum [MeV]",ylabel="Number of {0}".format(particle_full),plotname="momentum_{0}_{1}microns".format(particle,thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=None,ylow=None,yhigh=None,useTrueEdges=True)

    # Collect and save 2d histograms
    momentum_2d = infile.Get("momentum_xy_{0}".format(particle))
    momentum_2d.SetDirectory(0)
    thisThickness["momentum_xy_{0}".format(particle)] = momentum_2d
    position_2d = infile.Get("position_xy_{0}".format(particle))
    position_2d.SetDirectory(0)
    thisThickness["position_xy_{0}".format(particle)] = position_2d

    # Plot for visual fun
    myPainter.draw2DHist(momentum_2d,"momentum_xy_{0}_{1}".format(particle,thickness),"px [MeV]",-30,30,"py [MeV]",-30,30,particle_full,logz=True)
    myPainter.draw2DHist(position_2d,"position_xy_{0}_{1}".format(particle,thickness),"x [mm]",-1000,1000,"y [mm]",-1000,1000,particle_full,logz=True)

  savehists[thickness] = thisThickness

  # Make stacks of output particles within thickness
  to_stack = [thisThickness["polar_scattering_norm_neutron"],thisThickness["polar_scattering_norm_eplus"],thisThickness["polar_scattering_norm_gamma"],thisThickness["polar_scattering_norm_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [rad]",ylabel="Norm. particles per rad^{2}",plotname="stacked_particles_norm_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1,ylow=None,yhigh=None)

  to_stack = [thisThickness["polar_scattering_neutron"],thisThickness["polar_scattering_eplus"],thisThickness["polar_scattering_gamma"],thisThickness["polar_scattering_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [rad]",ylabel="Total particles",plotname="stacked_particles_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=1,ylow=None,yhigh=None)

  to_stack = [thisThickness["polar_scattering_deg_norm_neutron"],thisThickness["polar_scattering_deg_norm_eplus"],thisThickness["polar_scattering_deg_norm_gamma"],thisThickness["polar_scattering_deg_norm_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [deg]",ylabel="Norm. particles per deg^{2}",plotname="stacked_particles_deg_norm_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=90,ylow=None,yhigh=None)

  to_stack = [thisThickness["polar_scattering_deg_neutron"],thisThickness["polar_scattering_deg_eplus"],thisThickness["polar_scattering_deg_gamma"],thisThickness["polar_scattering_deg_eminus"]]
  myPainter.drawStackedHistos(to_stack,data=None,signal_list=None,stack_labels=["Neutrons","Positrons","Photons","Electrons"],data_label="",xlabel="Polar angle [deg]",ylabel="Total particles",plotname="stacked_particles_deg_{0}microns".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=90,ylow=None,yhigh=None)

# Close input file
infile.Close()

# Make plots comparing thicknesses
histos_list = [savehists[1]["polar_scattering_deg_eminus"],savehists[5]["polar_scattering_deg_eminus"],savehists[10]["polar_scattering_deg_eminus"]]
myPainter.drawOverlaidHistos(histos_list,data=None,signal_list=None,histos_labels=["1 #mum foil","5 #mum foil","10 #mum foil"],data_label="",xlabel="Polar angle [deg]",ylabel="Total electrons",plotname="electron_scattering_deg_compareThickness".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=45,ylow=None,yhigh=None)
histos_list = [savehists[1]["polar_scattering_deg_norm_eminus"],savehists[5]["polar_scattering_deg_norm_eminus"],savehists[10]["polar_scattering_deg_norm_eminus"]]
myPainter.drawOverlaidHistos(histos_list,data=None,signal_list=None,histos_labels=["1 #mum foil","5 #mum foil","10 #mum foil"],data_label="",xlabel="Polar angle [deg]",ylabel="Norm. electrons per deg^{2}",plotname="electron_scattering_deg_norm_compareThickness".format(thickness),doRatio=False,ratioName="",doBkgErr=False,logx=False,logy=True,nLegendColumns=1,extraLines=[],xlow=0,xhigh=45,ylow=None,yhigh=None)

# Save histograms that I created here, so I can look at them further if desired
outfile = ROOT.TFile.Open("plots.root","RECREATE")
outfile.cd()
for thickness in savehists.keys() :
  for name,hist in savehists[thickness].items() :
    hist.Write()
outfile.Close()


# Beam spread results, electrons only - compare methods
for thickness in 1, 5, 10 :
  # validate RMS calculation with one that is for sure OK
  test1 = savehists[thickness]["position_x_eminus"].GetRMS()
  test2 = compute_rms_1D(savehists[thickness]["position_x_eminus"])
  # print("Test RMS:",test1, test2) - validated
  test3 = savehists[thickness]["momentum_y_eminus"].GetRMS()
  test4 = compute_rms_1D(savehists[thickness]["momentum_y_eminus"])
  # print("Test RMS:",test3, test4) - validated

  # Not using ROOT's RMS since it compares to the mean along the x axis and I want to compare to zero
  # RMS of normalised distribution
  RMS_norm_dist = compute_rms_1D(savehists[thickness]["polar_scattering_norm_eminus"])
  print("RMS from solid-angle-normalised distribution {0:1.2g} = {1:1.2g} degrees".format(RMS_norm_dist,math.degrees(RMS_norm_dist)))
  # RMS by hand using bin center - non-normalised distribution
  RMS_notnormalised = compute_rms_1D(savehists[thickness]["polar_scattering_eminus"])
  print("RMS from unnormalised distribution {0:1.2g} = {1:1.2g} degrees".format(RMS_notnormalised,math.degrees(RMS_notnormalised)))

  # RMS of 2D distribution - position - and calculate from there
  RMS_position_2D = compute_rms_2D(savehists[thickness]["position_xy_eminus"])
  # Average using this is arctan RMS over 2 m
  RMS_angle_2D = math.atan(RMS_position_2D/2000.)
  print("RMS angle from 2D position distribution = {0:1.2g} = {1:1.2g} degrees".format(RMS_angle_2D,math.degrees(RMS_angle_2D)))

  # 1D (along each axis) of 2D distribution; how does that compare?
  RMS_xonly = math.atan(savehists[thickness]["position_xy_eminus"].GetRMS(1)/2000.)
  RMS_yonly = math.atan(savehists[thickness]["position_xy_eminus"].GetRMS(2)/2000.)
  print("For comparison, RMS along x only and y only give: {0:1.2g} {1:1.2g} = {2:1.2g} {3:1.2g} degrees".format(RMS_xonly,RMS_yonly,math.degrees(RMS_xonly),math.degrees(RMS_yonly)))