## Environment setup

On a laptop or similar, set up your Geant first:
```
cd /opt/local/libexec/Geant4/Geant4.10.6/
source geant4.sh
```
[go back to build dir]
```export G4LEDATA="/opt/local/share/Geant4/Data/Geant4.10.6/G4EMLOW7.9"```
[above corrects a mistake in the setup script in this release that will cause crashes]
```cmake -DGeant4_DIR=/opt/local/lib/Geant4/Geant4.10.6/Geant4-10.6.3 ..```

On Cedar, this setup script works:
```
module load StdEnv/2020
module load gcc/9.3.0
module load geant4/10.06
module load gnuplot
module load root

# Needed for CMake compilation
export G4DIR=/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Compiler/gcc9/geant4/10.06/lib64/Geant4-10.6.0
```

## Code setup

The following assumes you have defined `G4DIR` as above; otherwise replace with full path.
```
git clone https://github.com/kpachal/ScatteringRadiation.git
cd ScatteringRadiation
mkdir build; cd build
cmake -DCMAKE_CXX_STANDARD=17  -DGeant4_DIR=$G4DIR ..  
make -j4
```

This also creates a script in `build/` called `setup.sh`: you can source this to add the executables to your path so you can run them from any other location. I like this so I can use a run directory that isn't my build directory easily.

## Running

### Generating events

This is handled by `run_beam`. A specified number of electrons are created aiming towards a 2cm square tantalum target and are scattered, radiate, etc. At present (and later, by default) all electrons aim at exactly (0,0) so there is no beam spot size. This means outputs can be interpreted as an angular scattering probability per electron.

Relevant command line options are:
   * `nEvents` is number of beam electrons to generate (default 10000)
   * `targetThickness` in microns (default 1 micron)
   * `saveHighAngleOnly` to only save outgoing particles at least 20 degrees away from the beamline (default false)
   * `saveOnly` to only store specified particles in output file (default none). Currently supported options are `eplus`, `eminus`, `gamma`, or `neutron`
   * `seed` to set an explicit random number generator seed (useful for parallelisation).
   * `output` to set the string at the start of the output file name (will use "output" by default but add a few descriptors if non-default run options chosen)
   * `visuals` pulls up an interactive Geant window instead of following the other run directions. You can then generate some number of events and look at them using `/run/beamOn 100` or however many you want. At least for me, this doesn't then continue to an output file. Note: I am not sure if/how this works on Cedar as I have only used it locally.

Additional options exist but are not fully functional - for example adding more complex beam shapes.

An example run command to generate 1e5 electrons with a 2 micron foil, saving only photons at high angles, would be:
```
run_beam --nEvents 100000 --targetThickness 2 --saveHighAngleOnly --saveOnly gamma
```

## Analyzing outputs

The executable `analyze_data` can be run on any output file from the previous step to do a simple analysis in Root of the output and generate histograms. For example, to run on the output from the previous command,

```
analyze_data --input output_2micron_1e5events_30MeV_gamma.root
```

Additional command line options are:
   * `output` if you want the output file to be called something other than what's automatically generated
   * `tree` in case for some reason the input file no longer has tree name `Events` (you probably don't ever need this)