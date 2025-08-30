# pythia-wrapper
A complete icetray module that interfaces with PYTHIA8. It simulates a PYTHIA event and stores the particles into an IceCube IceTray I3MCTree

## Example set up file to correctly prepare CMakeLists
Note: pythia should be installed in your build directory for this example
```
#!/usr/bin/env bash
set -euo pipefail

# 1) Source the IceTray CVMFS environment
eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/setup.sh`
source /cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_7_x86_64_v2/bin/geant4.sh

# 2) Pythia/
export PYTHIA8_ROOT="$HOME/i3/icetray/pythia8315"
export PYTHIA8_INCLUDE_DIR="$PYTHIA8_ROOT/include"
export PYTHIA8_LIBRARY_DIR="$PYTHIA8_ROOT/lib"
export PYTHIA8_LIBRARY="$PYTHIA8_LIBRARY_DIR/libpythia8.so"
export PYTHIA8_LHAPDF6_LIBRARY="$PYTHIA8_LIBRARY_DIR/libpythia8lhapdf6.so"
export LD_LIBRARY_PATH="$PYTHIA8_LIBRARY_DIR:$LD_LIBRARY_PATH"


# 3) LHAPDF
export LHAPDF_LIBRARY_DIR=/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_7_x86_64_v2/spack/opt/spack/linux-centos7-x86_64_v2/gcc-13.3.0/lhapdf-6.5.1-6tgqvgixuyzyc2pgbel7krcl5nltte7i/lib
export LHAPDF_INCLUDE_DIR=/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_7_x86_64_v2/spack/opt/spack/linux-centos7-x86_64_v2/gcc-13.3.0/lhapdf-6.5.1-6tgqvgixuyzyc2pgbel7krcl5nltte7i/include

# 4) GEANT4
export GEANT4_DATA_DIR=/cvmfs/icecube.opensciencegrid.org/py3-v4.4.2/RHEL_7_x86_64_v2/spack/opt/spack/linux-centos7-x86_64_v2/gcc-13.3.0/geant4-data-11.2.2-ahaf4jla3qn7zpjde64hkp36tnd2a5ci/share/geant4-data-11.2.2

# 4) Re-cmake IceTray
~/i3/icetray/build/env-shell.sh
```

