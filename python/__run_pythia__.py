#!/usr/bin/env python 3

# Import libraries

import os, ctypes, sys, inspect, icecube
from os.path import expandvars, dirname, abspath, join, expanduser
from icecube import icetray, dataclasses, phys_services, dataio, 
from icecube.icetray import I3Logger

I3Logger.global_logger.set_level(I3LogLevel.LOG_INFO)
# ---------------------------------

load("pythia-wrapper", False)
del load

rng = phys_services.I3SPRNGRandomService(0, 1, 0)

def main():

    tray = I3Tray()

    cmd = os.path.expandvars("$I3_BUILD") + "/pythia-wrapper/resources/cmnd/neutrinoElectroweak.cmnd"

    print("Using cmnd file:", cmd, "exists?", os.path.exists(cmd))

    tray.AddModule("I3InfiniteSource", "src",
                  Stream=I3Frame.DAQ)

    tray.AddModule("pythia_wrapper", "pythia-wrapper",
                    PythiaCommandPath = cmd)

    # change name of pythia event output file here
    tray.AddModule("I3Writer", "hit_writer",
               Filename="pythia_output_10000_100TeV.i3")

    # change number of pythia events generated here
    tray.Execute(10000)
    print("Finished executing tray.")
    tray.Finish()

if __name__ == '__main__':
    main()
