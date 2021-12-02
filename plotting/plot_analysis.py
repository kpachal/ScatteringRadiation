import ROOT

# Prediction for electron destinations (multiple scattering): goudsmit-saunderson

infilename = "results_{0}micron_1e8events.root"

def normalise_solid_angle(inhist) :

    newhist = inhist.Clone(inhist.GetName()+"_norm")
    newhist.SetDirectory(0)
    newhist.Reset()

    # Solid angle subtended = 2 pi [-cos theta] | theta 1 to theta 2
    for bin in range(newhist.GetNbinsX())

for thickness in 1, 5, 10 :

  infile = ROOT.TFile.Open(infilename.format(thickness))
  polar_scattering_e = infile.Get("angle_polar_deg_eminus")
  polar_scattering_e.SetDirectory(0)

