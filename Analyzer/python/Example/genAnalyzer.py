# Load PyROOT
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gROOT.SetBatch(True)

# Load PyFWLite
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.FWLiteEnabler.enable()
from DataFormats.FWLite import Handle, Events

# other utilities
import sys

# define handles to read things from the event
gen_handle = Handle("std::vector<reco::GenParticle>")

fname = sys.argv[1] # take the name of the file to imput from the command line
events = Events(fname)
print "Reading %s " % fname,
for event in events:
    print "\n", ("--" * 75)
    print "Processing run %d, lumi %d, event %d" % (event.eventAuxiliary().run(), event.eventAuxiliary().luminosityBlock(), event.eventAuxiliary().event())

    # Get generated particles
    event.getByLabel("genParticles", gen_handle)
    # Make a list of those we are interested in (for a start, all within the detector coverage)
    gens = [ gp for gp in gen_handle.product() if abs(gp.eta()) < 5 ]
    # Print a numbered list of all them:
    print "Generated particles:"
    for igen, gen in enumerate(gens):
        print " %2d: gen   pt %7.2f  eta %+5.2f  phi %+5.2f  energy %7.2f   pdgId % +5d  charge %+1d " % (igen, gen.pt(), gen.eta(), gen.phi(), gen.energy(), gen.pdgId(), gen.charge())
    print ""
