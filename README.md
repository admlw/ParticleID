# Particle ID module

This module is intended to be used for analyses on MicroBooNE. It's currently under active development.

The Analyzer provided in the Modules directory has one dependency, the UBXSec module, which can be found here: https://github.com/marcodeltutto/UBXSec. This allows running over events selected by the CC-Inclusive analysis to get selected tracks, which are more likely to be protons than just running over all tracks in the TPC. This dependency can be removed by removing the relevant include statements and the relevant lines in the Modules/CMakeLists.txt file.
