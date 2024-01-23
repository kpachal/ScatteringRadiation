import ROOT
import math
from art.morisot_2p0 import Morisot_2p0

# Initialize painter
myPainter = Morisot_2p0()
myPainter.colourpalette.setColourPalette("notSynthwave")
myPainter.labelType = 4
myPainter.savePDF = True
myPainter.luminosity = -1
myPainter.CME = -1

infile = "../results_1micron_1e7events_30MeV.root"

openfile = ROOT.TFile.Open(infile,"READ")

hist_x = openfile.Get("angle_x_eminus")
hist_x.SetDirectory(0)
hist_y = openfile.Get("angle_y_eminus")
hist_y.SetDirectory(0)

# Now fit each
for axis,histo in zip(["x","y"],[hist_x,hist_y]) :
    gfit = ROOT.TF1("fit{0}".format(histo.GetName()),"gaus",-0.1,0.1)
    histo.Fit(gfit,"L","",-0.02,0.02)
    pars = gfit.GetParameters()
    print(pars[0],pars[1],pars[2])

    # Plot with fit over full range
    myPainter.drawHistsWithTF1s([histo],[gfit],hist_labels=["g4 sim distribution"],func_labels=["Gaussian fit"],xlabel="Angle in {0} [rad]".format(axis),ylabel="Electrons",plotname="plots/fitted_{0}_angle".format(axis),doRatios=False,ratioName="Blah",nLegendColumns=1,extraLines={"Fit range between -0.02 and 0.02"},xlow=-0.1,xhigh=0.1)