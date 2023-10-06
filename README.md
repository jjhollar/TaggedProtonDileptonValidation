# TaggedProtonDileptonValidation

****** Introduction ******
Code for performing a single-arm dimuon selection with PPS protons reconstructed 
on-the-fly, using new alignment+optics settings. 

This is a heavily cleaned up and updated version of the code used for Run 2. 
So far it is only tested to compile & run in CMSSW_13_0_10, on Run 2 (2018) data. 

****** Overview ******
The code consists of 2 parts: 
 1) A CMSSW analyzer that runs on AOD, with the option to re-reconstruct the protons. 
    The output is a plain ROOT ntuple with many branches containing information on 
    the muons, protons, PPS local tracks, and vertex/track counting
 2) A ROOT macro that runs on the output of step 1), applies some "analysis" cuts, 
    and outputs final distributions. The output includes a ROOT file with histograms, 
    and a text file with information on the signal-like events

****** Setup and running the analyzer ******

       cmsrel CMSSW_13_0_10

       cd CMSSW_13_0_10/src/

       cmsenv

       git clone https://github.com/jjhollar/TaggedProtonDileptonValidation.git
       (note this will be moved to a more central location in the POG repository)

       scram b

       cd TaggedProtonDileptonValidation/GammaGammaLeptonLepton/test/       

       cmsRun RunGammaGammaLeptonLepton_cfg.py       

The config file RunGammaGammaLeptonLepton_cfg.py should be modified to choose the 
input files, number of events, trigger names, etc. The variable "DORERECO" controls 
whether the analyzer uses pre-existing protons from the AOD central production, or 
runs the proton reconstruction on-the-fly and uses the new objects. 

The global tag and conditions should be set as needed to pick up the new alignment/optics 
setting, depending on the format in which they're provided (in a tag, or sqlite file, or other)

****** Running the ROOT macros ******

An example ROOT macro to make final histograms is in TaggedProtonDileptonValidation/GammaGammaLeptonLepton/macros
The file DimuonMacro.h should be modified to use the output files from the analyzer as input.

    root -l 
    root [0] .L DimuonMacro.C
    root [1] DimuonMacro m
    root [2] m.Loop()

The Loop function takes several arguments that can be used to configure the job:
    multi: use events with a single proton per arm (=0) or allow multuple protons (=1)
    mc: include some generator-level information when running on mc (=1) instead of data (=0)
    sb: apply sideband-region cuts in acoplanarity and Ntracks (=1) instead of signal-region cuts (=0)
    yr: year of data analyzed. This is less important in Run 3, where only pixel tracking is used
    nearfar: for single-RP distributions, choose between the near/210m pixels (=1) or far/220m pixels (=0)

The cuts applied to the muons, and track multiplicity at the dimuon vertex, are copied from the original 2016 
analysis, and may not be optimal for Run 3. 