
BEAST 2 (http://beast2.org) package for logging innovations on the branch leading to a given node.

The Java package in the innovations BEAST2 package is:

### `innovations.evolution.likelihood`
* 'InnovationLogger' - takes as input a taxonset, then performs an ancestral state reconstruction. For a given character, if the parent node of the MRCA of the taxonset is reconstructed as having an absent state, and the MRCA node as having a present state, the logger will record 1 for that character. Otherwise 0 will be recorded. The columns of the log file then just need to be averaged to get the posterior probability of an innovation in each character for the taxonset.

There is an R function in the [beast_asr](https://github.com/king-ben/beast_asr) repo for automating the generation of the relevant xml code for this logger.

See `examples/make_innovation_blocks.R` for an example script for making the xml input, including an explanation of all the options.

See `examples/analyse_innovations.R` for a script that analyses the output of the innovations log file.

See `examples/Malekula_asr.xml` for a simple example beast2 xml file.

Installing package manually
---------------------------
Ensure you are running at least beast 2.7.
Download the zip file from [releases](https://github.com/king-ben/innovations/releases). Copy this zip file to the BEAST addon directory in your computer and unzip it. To find out where the addon directory is, open beauti, File - manage packages, click on the question mark in the bottom right corner. You may also need to clear the beast class path. This is possible in the file menu of beauti.
