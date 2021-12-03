import ROOT
import math

# Prediction for electron destinations (multiple scattering): goudsmit-saunderson

infilename = "../results_{0}micron_1e8events.root"

def normalise_solid_angle(inhist) :

    newhist = inhist.Clone(inhist.GetName()+"_norm")
    newhist.SetDirectory(0)
    newhist.Reset()

    # Solid angle subtended = 2 pi [-cos theta] | theta 1 to theta 2
    for bin in range(newhist.GetNbinsX()+1) :
        fullval = inhist.GetBinContent(bin)
        theta1 = inhist.GetBinLowEdge(bin)
        theta2 = inhist.GetBinLowEdge(bin+1)
        normval = 2 * math.pi * (- math.cos(theta2) + math.cos(theta1))
        newhist.SetBinContent(bin,fullval/normval)
        newhist.SetBinError(bin,inhist.GetBinError(bin)/normval)

    return newhist

savehists = []
for thickness in 1, 5, 10 :

  thisname = infilename.format(thickness)
  print("Opening file",thisname)
  infile = ROOT.TFile.Open(thisname)
  polar_scattering_e = infile.Get("angle_polar_eminus")
  polar_scattering_e.SetDirectory(0)
  savehists.append(polar_scattering_e)

  polar_scattering_e_norm = normalise_solid_angle(polar_scattering_e)
  polar_scattering_e_norm.SetName(polar_scattering_e_norm.GetName()+"_{0}micron".format(thickness))
  savehists.append(polar_scattering_e_norm)
  
outfile = ROOT.TFile.Open("plots.root","RECREATE")
outfile.cd()
for hist in savehists :
    hist.Write()
outfile.Close()