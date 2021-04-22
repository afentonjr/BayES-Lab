# BayES’Lab (Bayesian Electroactive Species Labeling)

BayES’Lab (Bayesian Electroactive Species Labeling) is a repository of functions that can both 1) estimate electrochemical-transport parameters of redox active species and 2) infer the identities of electroactive compounds catalogued in a library from an experimental electrochemical dataset. Each of these functions is separately presented as its own subroutine and can be run independently, if desired.

See below for implementation instructions and details. Further, this repository was developed for and in tandem with the manuscript “Using voltammetry augmented with physics-based modeling and Bayesian hypothesis testing to estimate electrolyte composition" (to be submitted). <b><i>Please cite this paper if you use this package.</b></i>

<b>Code style:</b> MATLAB®

## Overview
  
Voltammetric capabilities may be enhanced by joining inferential analysis with electrochemical physical modeling. One such application is automated identification of electroactive compounds. As such, the goal of this project was to use first-principles simulation and Bayesian inference to label electroactive compounds by 1) creating a library with training datasets and 2) applying this library to testing data.

This labeling process is achieved with two corresponding and overarching modules: <b>Library Development</b> and <b>Compound Identification</b>. As already mentioned, both are included in this repository and can be run either sequentially (<b>Library Development</b> --> <b>Compound Identification</b>) or separately.

The <b>Library Development</b> module takes multiple cyclic square wave (CSW) voltammograms and other relevant inputs (e.g., experimental parameters, initial guesses) to estimate the electron transfer mechanism and the corresponding electrochemical-transport descriptors. <b>Although the intended purpose of this code was to create a library of species for compound identification, it can also be used to extract electrochemical-transport descriptors from a known compound of interest for its own end</b> using CSW voltammograms acquired on a planar working electrode in a three-electrode configuration.

The <b>Compound Identification</b> module, in turn, uses a library of parameters (either externally supplied or assembled using the <b>Library Development</b> module) to examine new experimental data and estimate the composition of redox-active compounds in the electrolyte. <b>The Compound Identification module can be used to identify electroactive species catalogued in a library using either <i>in situ</i> CSW voltammetry or cyclic voltammetry (CV) in near-real time (> 3 min).</b>

## Framework

This code runs on MATLAB® R2020a with the Parallel Computing Toolbox optionally installed.

## Installation

Information on how to install MATLAB® may be found <a href="https://www.mathworks.com/help/install/install-products.html">here</a>.
The Parallel Computing Toolbox, in turn, can be found <a href="https://www.mathworks.com/products/parallel-computing.html">here</a>.

## Use

The following sections contain information on how to use BayES’Lab. As mentioned above, there are two modules: <b>Library Development</b> and <b>Compound Identification</b>. These may be run either indepdently or sequentially; for the latter, a library can be constructed using the <b>Library Development</b> module, which can then be applied to the <b>Compound Identification</b> module to identify electroactive compounds. Note that for each module, there are a list of inputs/outputs, example script(s), command-line use instructions, and for <b>Library Development</b>, instructions for use on the MIT Supercloud (in theory, any supercomputing resource can be used with the appropriate modifications).

### Library Development

#### General

This module estimates both the electron transfer mechanism, along with the corresponding electrochemical-transport descriptors, from experimental CSW voltammograms.

<b>Inputs</b>
- <b>data:</b> compiled CSW voltammograms from which parameters and the electron transfer mechanism will be estimated
- <b>paramsin:</b> input parameters necessary for the module to operate. Examples include the working electrode radius, the temperature, and the solution resistance
- <b>paramguess:</b> initial guesses for the electrochemical-transport descriptors to be estimated
- <b>conc:</b> concentrations of the species that comprise each CSW voltammogram (mol L<sup>-1</sup>)
- <b>ninitguess:</b> number of initial guesses to be conducted
- <b>bsubtract:</b> the difference current when no redox-active species are present; used for background correction (&mu;A)
- <b>data_ref:</b> CSW voltammogram of the redox-active compound of interest and a second compound used to calibrate a pseudoreference electrode to a known redox event
- <b>data_bgref:</b> CSW voltammogram of only the compound used to calibrate the pseudoreference; used to calibrate the background current
- <b>prior:</b> prior probabilities for each electron transfer mechanism. Used in Bayes' Rule

<b>Outputs</b>
- <b>paramsfound:</b> the estimated parameters
- <b>paramslb:</b> the lower bound of the estimated parameters
- <b>paramsub:</b> the upper bound of the estimated parameters
- <b>ETM:</b> the estimated electron transfer mechanism (1 = reversible, 2 = quasireversible kinetics)
- <b>PETM:</b> the probability of the selected electron transfer mechanism
- <b>dataEstep:</b> the post-corrected potential data, averaged across all input CSW voltammograms (V vs. ref)
- <b>dataInet:</b> the post-corrected difference current, averaged across all input CSW voltammograms (A)
- <b>inetsim:</b> the best-fit simulated current for all electron transfer mechanisms; one column per ETM (A)
- <b>Eref:</b> the redox potential of the calibration redox event with respect to the pseudoreference. Currently reports for the final voltammogram in the input dataset (V vs. ref)
- <b>simtime:</b> simulation time (s)

#### Example Script

To run an example script for the library development of 1 mM unsubstituted phenothiazine in 0.1 M TBAPF6 in dichloromethane, <b>download the repository, set the directory to be "BayES-Lab-test-v01/", open the script "Example_Script_Library_Development_PT.m", and run the script (no function inputs are needed).</b> Note that the progress of the optimizer is output onto the command window to display its progress.

#### Local Command Line Use (e.g., on a laptop)

Instructions will be provided to run the code using the data provided for unsubstituted phenothiazine (PT). Load the files for the respective folder by opening the files directly within the folder "Library Development/PT Inputs and Previous Results/". Then, load all the files in the PT folder that end with a ".mat" label. Conversely, you can navigate to the directory in the MATLAB® application and type “load(‘___.mat’)” into the command window for all the files.

We note that the automated data processing has only been verified for cases where the initial potential sweep direction is oxidative (i.e., increasing in working eletrode potential); performance when the initial sweep is reductive has not been verified. However, the parameter fitting routine is expected to work when the initial potential sweep is reductive if the data is manually processed.
 
Subsequently, exit the "PT Inputs and Prevous Results/" directory (i.e., the directory should now be "Library Development/") and type the following into the command window:

```

params_inPT{6}=200;
[paramsfoundPT,paramslbPT,paramsubPT,ETMPT,PETMPT,dataEstepPT,dataInetPT,inetsimPT,ErefPT,simtimePT]=library_entry_params_compound(dataPT,params_inPT,params_guessPT,conc_inPT,0,bsubtract,dataPTref,dataFc_trunc,priorPT);

```

Note that the time mesh is expanded to 5 ms (from 1 ms), and only 1 initial guess (=ninitguess+1) is used, to achieve a simulation time that is feasible on a local computing resource (<i>ca</i>. 1 h using two cores with minimal background processes). Also note that for different compounds (e.g., MPT, EPT), appropriate variable names will have to be modified within the input string (e.g., “dataPTref” must become “dataMPTref”, "params_inPT" must be "params_inMPT").

#### Supercomputer use

It may be desirable to run the "Library Development" protocol with a greater time resolution, more initial guesses, or both. To this end, it is also possible to run the "Library Development" module on a supercomputing resource. Within the published work, we performed parameter estimation using the MIT Supercloud supercomputing resource. <i>The following instructions have only been verified for the MIT Supercloud supercomputing resource; alterations in procedure may be necessary for other supercomputing resources</i>. To perform the routine, we used an emacs script shown below as an example for PT:

```

#!/bin/bash

#SBATCH -n 100
#SBATCH -N 4

matlab -nodisplay -singleCompThread -r 

"[pfoundPT,plbPT,pubPT,ETMPT,PETMPT,dataEstepPT,dataInetPT,inetsimPT,ErefPT,simtimePT]=library_entry_params_compound(load('2020-06-30-PT-Data-for-Param-Estimation.mat'),load('2020-06-30-Input-Params-for-PT.mat'),load('2020-06-30-Guess-Params-for-PT.mat'),load('2020-06-30-Input-Concentrations-for-PT.mat'),99,load('2020-06-30-Background-Data-for-Background-Subtraction.mat'),load('2020-06-30-PT-Data-with-Ferrocene-for-Param-Estimation-Potential-Reference.mat'),load('2020-06-30-Background-Data-for-Background-Subtraction-Potential-Reference.mat'),load('2020-07-21-Priors-for-PT-Parameter-and-ETM-Estimation.mat')); exit"

```

Note that each of the variables must be loaded and that 100 initial guesses are used (=ninitguess+1). The job can then be submitted to the HPC resource using the above emacs script. We note that with the given settings, the total simulation time for each phenothiazine is <i>ca</i>. 10 days, although the simulation time can be reduced by appropriate parallelization and by adjusting the time resolution or number of initial guesses.

Additional modifications must also be made to the MATLAB® code sent to MIT Supercloud. Specifically, the blocks of code below must be uncommented in order for the code to work on MIT Supercloud.


```
% bsubtract=bsubtract.bsubtract;
% conc=conc.conc_inPT;
% data=data.dataPT;
% data_bgref=data_bgref.dataFc_trunc;
% data_ref=data_ref.dataPTref;
% paramguess=paramguess.params_guessPT;
% paramsin=paramsin.params_inPT;
% prior=prior.priorPT;

% save('2021-02-02-Estimated-PT-Params-100ig-1000pps.mat','paramsfound')
% save('2021-02-02-Estimated-PT-Params-LB-100ig-1000pps.mat','paramslb')
% save('2021-02-02-Estimated-PT-Params-UB-100ig-1000pps.mat','paramsub')
% save('2021-02-02-Estimated-PT-ETM-100ig-1000pps.mat','ETM')
% save('2021-02-02-Estimated-PT-ETM-Prob-100ig-1000pps.mat','PETM')
% save('2021-02-02-Estimated-PT-Potential-Data-100ig-1000pps.mat','dataEstep')
% save('2021-02-02-Estimated-PT-Current-Data-100ig-1000pps.mat','dataInet')
% save('2021-02-02-Estimated-PT-Current-Simulated-100ig-1000pps.mat','inetsim')
% save('2021-02-02-Estimated-PT-Reference-Potential-100ig-1000pps.mat','Eref')
 
%% Stop timer
simtime=toc;
% save('2021-02-02-Estimated-PT-Sim-Time-100ig-1000pps.mat','simtime')

```

The first block of lines must be uncommented to convert structures to a usable format (vectors and matrices); files directly loaded into the function in the emacs script above are input as structures. The other blocks are uncommented to save the function outputs; MIT Supercloud does not save the function outputs otherwise. Note that the names of the variables will also have to be adjusted for different compounds (e.g., conc=conc.conc_inMPT for MPT, etc.).

These blocks correspond to lines 183-190, lines 230-238, and line 242, respectively, in library_entry_params_compound.m.

### Compound Identification

#### General

This module estimates the composition of redox-active species in an electrolyte, and currently accepts both CVs and CSWVs as inputs.

<b>Inputs</b>
- <b>data:</b> a single voltammogram from which compounds will be estimated
- <b>params:</b> input parameters necessary for the module to operate. Examples include the working electrode radius, the temperature, and the solution resistance
- <b>ninitguess:</b> number of initial guesses to be conducted
- <b>lib:</b> the library to be referenced for compound identification
- <b>bsubtract:</b> the (difference) current when no redox-active species is present; used for background correction (&mu;A for CSW voltammetry, mA for CV)
- <b>data_ref:</b> voltammogram of the redox-active compound of interest and a second compound used to calibrate a pseudoreference electrode to a known redox event
- <b>data_bgref:</b> voltammogram of only the compound used to calibrate the pseudoreference; used to calibrate the background current
- <b>prior:</b> prior probabilities for each electron transfer mechanism. Used in Bayes' Rule
- <b>type:</b> numerical input to determine whether a CV or CSW voltammetry is being analyzed

<b>Outputs</b>
- <b>P:</b> probability of each compound being present in the electrolyte
- <b>conc0:</b> the estimated concentrations of the compounds estimated using regression alone (mol m<sup>-3</sup>)
- <b>concfinal:</b> the estimated concentrations of the compounds estimated using both regression and binary hypothesis testing (should be more accurate than conc0, in units of mol m<sup>-3</sup>)
- <b>conclbfinal:</b> lower bound estimate of concfinal (mol m<sup>-3</sup>)
- <b>concubfinal:</b> upper bound estimate of concfinal (mol m<sup>-3</sup>)
- <b>Estep:</b> the post-corrected potential data (V vs. ref)
- <b>inet0:</b> the best-fitting (difference) current using regression alone (A)
- <b>inet:</b> the best-fitting (difference) current using both regression and binary hypothesis testing (A)
- <b>dataInet:</b> the post-corrected (difference) current (A)
- <b>Eref:</b> the redox potential of the calibration redox event with respect to the pseudoreference. Outputs 0 if no correction is desired (V vs. ref)
- <b>simtime:</b> simulation time (s)

#### Example Script

To run an example script to estimate the composition of an electrolyte containing unsubstituted and methyl phenothiazine using CV, <b>download the repository, open the directory to be "BayES-Lab-test-v01/", open the script "Example_Script_Compound_Identification_CV.m", and run the script (no function inputs are needed).</b> Note that similarly to the library development module, text is output onto the command window to display the progress of the protocol.

<b>To run an example script for CSW voltammetry taken from the same electrolyte, perform the above procedure for "Example_Script_Compound_Identification_CSWV.m".</b>

Note that with the exact inputs in this script, two warnings will appear that say there will be no adjustment of the reference potential. This is normal and expected. If the Parallel Computing Toolbox is being used, a warning may be issued for every type of parameter in the library, as they are uninitialized temporary variables; this is also okay and expected.

#### Local Command Line Use (e.g., on a laptop)

As with the library construction section, load the files by opening all the files with a ".mat" label in the "Compound Identification/" folder. Conversely, you can navigate to this directory in MATLAB® and type “load(‘___.mat’)” into the command window.
 

Subsequently, type the following input into the command line.

```

[P,conc0,concfinal,conclbfinal,concubfinal,Estep,inet0,inet,dataInet,Eref,simtime]=conc_compound_find(PTMPTdatafirstelyte(:,1:4),params_indata,99,library,bsubtract,0,dataFc_trunc,prior,0);

```

## Developer Use

Navigate to the directory where you would like the local copy of the source code to exist, and then clone the repository using:

```

git clone git@github.mit.edu:afenton/BayES-Lab.git

```

Conversely, you may directly download the repository from the online page.

Congratulations! BayES'Lab is now installed! With this basic installation, you will have access to the functions therein.

# Credits
Alexis M Fenton Jr. | Fikile R. Brushett | MIT Supercloud

# Contribute
Thank you for your interest in augmenting the capabilities of voltammetry through inferential methods!

The following guidelines can help you get started within this project. Use your best judgement and coding knowledge to help improve it, and please feel free to propose changes in a pull request.

## Getting started

BayES’Lab is run exclusively in MATLAB®. Refer to the Installation section to access and download MATLAB®.

## Issues

Before submitting your own issue, make sure the same one does not already exist by searching under <a href="https://github.mit.edu/afenton/BayES-Lab-test-v01/issues">issues</a>. If you do not locate a similar issue already, feel free to open a <a href="https://github.mit.edu/afenton/BayES-Lab-test-v01/issues/new">new issue</a>. When writing your issue, consider the following points and be as detailed as possible:

1.	Write step by step directions for replicating your issue
2.	Describe what you would expect to happen and what you actually get
3.  Provide relevant screenshots
4.	Include your operating system version and MATLAB® version

## Pull Requests

Pull requests are always welcome for suggestions to improve either the code or usability. Before submitting the pull request, please ensure that your standalone code is working properly by both running the existing tests and adding tests of any new functionality. Please explain any new tests and functionality in detail.

# API Reference
A general description of variables is provided here. More specificity, including the corresponding units to each variable, may be found in the header of each function.

## Functions

### Library Development

<b>library_entry_params_compound</b>(data,paramsin,paramsguess,conc,ninitguess,bsubtract,data_ref,data_bgref,prior)
- data = experimental data for which the electron transfer mechanism and the corresponding parameters are estimated. Must be a set (n>1) of CSW voltammograms
- paramsin = input parameters, such as the CSW voltammogram waveform and the temperature
- paramsguess = a guess for the parameters that will be estimated
- conc = a vector of concentrations, used for fitting the experimental data with a model
- ninitguess = number of initial guesses. The parameter space is assumed to be very noisy, with many local minima. As such, many initial guesses are likely needed to increase the changes of obtaining the global minimum 
- bsubtract = a scan with no redox couples present. Used to find the background current for subtraction from the experimental data
- data_ref = the experimental data with a reference redox couple added to adjust the potential to be referenced versus a desired redox event (e.g., ferrocene)
- data_bgref = a dataset where the only active species is the reference redox couple. Used to adjust the potential of the background current
- prior = prior probabilities of a given electron transfer mechanism
- return = paramsfound (estimated parameters), paramslb (lower bound of the estimates), paramsub (upper bound of the estimates), ETM (the estimated electron transfer mechanism), PETM (the probability of the estimated electron transfer mechanism), dataEstep (the processed experimental potential), dataInet (the processed experimental current), inetsim (the simulated current), Eref (the redox potential of the reference redox couple vs. the original reference electrode used), and simtime (the simulation time, for internal consumption)

Accepts experimental data and outputs the estimated electron transfer mechanism and the related parameters.


<b>data_processing_library_find</b>(data,bsubtract,data_ref,data_bgref,iter,option,paramsin,conc)
- data = same as in library_entry_params_compound
- bsubtract = same as in library_entry_params_compound
- data_ref = same as in library_entry_params_compound
- data_bgref = same as in library_entry_params_compound
- iter = indicator for the number of trials to be averaged in the parameter / ETM estimation process
- option = option to allow / not allow for background correction
- paramsin = same as in library_entry_params_compound
- conc = same as in library_entry_params_compound
- return = dataV (the processed experimental potential), dataInet (the processed difference current), Idatastd (the experimental standard deviation of the current), dataEstep (the baseline potential), and Eref (the amount the potential was adjusted based on that of the reference redox couple)

Pre-processes the data by performing iR correction, adjusting the potentials to that of a different, desired redox couple (if desired), and performing background subtraction. Enables data pre-processing in a consistent manner.


<b>iR_correction</b>(dataV,dataI,R,comp)
- dataV = the processed experimental potential
- dataI = the processed experimental current
- R = the solution resistance
- comp = the compensation automatically applied by the potentiostat
- return = dataEstep (the shifted baseline potential), dataVnew (the shifted potential)

iR-corrects potentials (both baseline and regular).


<b>reference_correction</b>(dataEstep_comp,dataInet_comp,dataV,data_ref,Rref,comp,bg)
- datatEstep_comp = the baseline potential of the voltammogram with the compound of interest
- dataInet_comp = the difference current of the voltammogram with the compound of interest
- dataV = the potential of the voltammogram with the compound of interest
- data_ref = same as in library_entry_params_compound
- Rref = the resistance of the electrolyte containing both the compound of interest and the reference redox couple
- comp = the compensation automatically implemented by the potentiostat
- bg = indicator of whether the background voltammogram or that with the compound of interest is having its potential shifted
- return = dataEstep (the shifted baseline potential), dataVnew (the shifted potential), Eref (the amount the potentials were shifted by)

This function shifts the potentials of the voltammograms such that it is referenced to a redox event of interest. For example, potentials may be shifted such that the redox event is 0 V vs. ferrocene.


<b>voltammogram_params_fitting_protocol</b>(dataV,dataInet,dataIstd,ninitguess,paramsin,paramguess,conc)
- dataV = same as in library_entry_params_compound
- dataInet = same as in library_entry_params_compound
- dataIstd = the experimental standard deviation of the current, calculated during automated data pre-processing
- ninitguess = same as in library_entry_params_compound
- paramsin = same as in library_entry_params_compound
- paramguess = same as in library_entry_params_compound
- conc = same as in library_entry_params_compound
- return = params (the estimated parameters), inetsim (the simulated difference current using the optimal set of parameters), lowerbound (the lower bound of the parameters), and upperbound (the upper bound of the parameters)

This function inputs processed experimental data to estimate an electron transfer mechanism and the relevant parameters.


<b>biccompx</b>(dataInet,dataIstd,inetsim,data)
- dataInet = processed experimental difference current
- dataIstd = experimental standard deviation of the difference current
- inetsim = the simulated difference current for the optimal set of parameters
- data = same as in library_entry_params_compound
- return = BIC, the value of the Bayesian Information Criterion for a given electron transfer model

This function calculates the BIC for a given electron transfer model, used for comparison with others to determine the most likely model.


<b>obj_func_min_param_fit_CSWV</b>(c,params_fit,dataInet,model,dataV,dataIstd)
- c = the vector of parameters to be estimated
- params_fit = fixed parameters used to simulate the model
- dataInet = the processed experimental difference current
- model = the electron transfer mechanism being considered (e.g., reversible, quasireversible)
- dataV = the processed experimental potential
- dataIstd = the experimental standard deviation of the difference current
- return = f (the objective function value to be minimized)

This function finds the parameters that maximizes the likelihood function. It is assumed that the likelihood function is a multivariate normal distribution, which means that this maximization is performed by minimizing the argument in the exponent of the likelihood function.


<b>obj_func_find_param_bounds</b>(c,params_fit,dataInet,model,dataV,dataIstd,fval,params,q)
- c = the vector of parameters to be estimated
- params_fit = fixed parameters used to simulate the model
- dataInet = the processed experimental difference current
- model = the electron transfer mechanism being considered (e.g., reversible, quasireversible)
- dataV = the processed experimental potential
- dataIstd = the experimental standard deviation of the difference current
- fval = the minimum value of the exponential argument of the multivariate normal distribution found in obj_func_min_param_fit_CSWV
- params = the best estimate of the parameters
- q = an index used to determine which parameter the bounds for which the bounds are being estimated
- return = f (the objective function value to be minimized)

This function is used to estimate the lower and upper bounds of the parameters by finding the value of the parameter of interest which halves the likelihood function (an empirical standard) when all other parameters are kept at their optimal values.


<b>Voltammetry_Simulator_MacroDisk_2019_11_07</b>(kinetics,direction,type,dataV,params)
- kinetics = indicates whether reversible or quasireversible electron transfers are being considered
- direction = an indicator of whether the initial voltammetric sweep is oxidative or reductive
- type = indicates whether a conventional cyclic voltammogram (CV) or CSW voltammogram is being simulated
- dataV = a vector of processed experimental potentials
- params = relevant input parameters
- return = conc (output concentrations as a function of space and time), potential (the output potential), I (the output current), t (the output time vector), phiLSE (the concentration normalized current, for the compound identification protocol)

This function simulates either a conventional CV or a CSW voltammogram using a finite difference numerical framework.


<b>E_find_interpolate</b>(dataV,sub,tau)
- dataV = the processed experimental potential
- sub = the number of desired steps in the simulated voltmmogram (e.g., if there are 1000 time steps, sub=1000)
- tau = the time duration of each potential hold
- return = t (the expanded time vector), E (the expanded potential vector)

This function interpolates the potential to the desired number of substeps. This is useful for simulating voltammograms with small time steps


<b>solveconc</b>(X,params,theta,cinit,h)
- X = the dimensionless spatial mesh
- params = imported parameters
- theta = dimensionless overpotential
- cinit = the initial concentration guess for the concentration vector to be solved. If a voltammogram is being simulated across time, this may be the concentration vector at the previous time point
- h = the initial dimensionless step size
- return = conc (the concentration profile of the active species at a given time point during the simulation of a voltammogram)

This function sets up and solves the governing partial differential equations in a voltammetric simulation by creating a finite difference framework and solving via linear algebra.


<b>Difference_current_extract</b>(J,E,b,ind,maxind)
- J = the dimensionless flux to the surface of the electrode
- E = a vector with the dimensional potential. The first entry is the dimensional potential at the previous time point, and the second is the dimensional potential at the current time point
- b = the time index for the beginning of the current pulse, not including the time point about to be evaluated. Used for tracking the pulse duration
- ind = the current time index
- maxind = the maximum index possible
- return = Jdiff (a vector of dimensionless fluxes that are marked to either be kept or deleted in the final calculation of the flux / current), Ediff (the analogous case of Jdiff for potential), b (the time index for the beginning of the current pulse, including the time point evaluated in the function)

This function serves as the first of two functions to increase the time step size of the simulated voltammogram to be the same as the input experimental data. The overall process achieves this by selectively choosing which output values not to report. This function marks at which time points the dimensionless fluxes and dimensional potentials should be eliminated.


<b>difference_current_calculate</b>(Ediff,Jdiff)
- Ediff = the vector of dimensional potentials marked for keeping or deletion by Difference_current_extract
- Jdiff = the vector of dimensionless fluxes marked for keeping or deletion by Difference_current_extract
- return = Estep (the compressed vector of dimensional potentials), Jnet (the compressed vector of dimensionless fluxes)

This function compresses the vectors of potentials and fluxes by deleting the time points marked by Difference_current_extract




### Compound identification


We note that multiple functions in this module are similar, but not the same, as that in the Library Development module (e.g., data pre-processing). We choose not to include these functions here for brevity, rather including only new parent functions. For details on all functions, download or clone the code for these parent functions.


<b>conc_compound_find</b>(data,params,ninitguess,lib,bsubtract,data_ref,data_bgref,prior,type)
- data = experimental data for a single voltammogram (either a CV or a CSW voltammogram)
- params = relevant parameters for the experiment, such as the radius of the working electrode and the voltammetric waveform (e.g., pulse height, scan rate)
- ninitguess = the number of randomized initial concentration guesses the protocol will iterate through
- lib = library of compounds, with relevant parameters listed
- bsubtract = a dataset of blank electrolyte (only solvent and supporting salt) for baseline subtraction. Multiple trials are necessary to evaluate a mean background current and its standard deviation
- data_ref = a single experimental voltammogram of the electrolyte of interest that also has added a redox-active compound whose redox potential can be used for referencing the applied potential. An example is an electrolyte with ferrocene intentionally added to reference all the redox events to the redox potential of ferrocene
- data_bgref = multiple trials of a blank electrolyte containing a redox-active compound used for referencing the potential. Used both to reference the potential of the background current and to estimate the standard deviation of the current as a function of the magnitude of the current
- prior = prior probability of existence assigned to each compound in the library, in line with the Bayesian framework
- type = determines whether cyclic voltammetry data or CSW voltammetry data is being analyzed. type=1 means a CV is being analyzed; any other number means a CSW voltammogram is being analyzed
- return = P (the probability that each compound is present in the electrolyte), conc0 (the vector of best-fit concentrations for all the compounds in the library before culling via Bayesian inference), concfinal (the vector of best-fit concentrations after culling occurs), conclbfinal (the lower bound estimate of concfinal), concubfinal (the upper bound estimate of concfinal), Estep (the post-processed experimental potential or baseline potential), inet0 (the vector of simulated currents before culling), inet (the vector of simulated currents after culling), dataInet (the post-processed data current to be fitted), Eref (the redox potential of the reference redox couple vs. the original reference electrode used), simtime (the time for the simulation)

This function is the parent function used to estimate the probabilities of each compound being present in solution. It first pre-processes the data (e.g., ohmic compensation), and then simulates concentration-normalized voltammograms based on the library parameters for each compound. These concentration-normalized voltammograms are then regressed to a new experimental dataset to yield a vector of best-fit concentrations. Identification of each compound then takes place by accounting for the goodness of fit and the parsimoniousness of the model when each compound is excluded (the null hypothesis) and when it is included (the alternative hypothesis) to report the probability that the compound in question is present.


<b>voltammogram_conc_fitting_protocol</b>(dataInet,dataIstd,phi,ninitguess,concguess,lib,fittype)
- dataInet = post-processed data current
- dataIstd = the estimated standard deviation in the current for the experimental data of interest
- phi = the simulated concentration-normalized current
- ninitguess = same as in conc_compound_find
- concguess = an initial guess for the concentrations
- lib = same as in conc_compound_find
- return = conctemp (vector of best-fit concentrations of all compounds in the library), lbfind (the vector of estimated lower-bound best-fit concentrations), ubfind (the vector of upper-bound best-fit concentrations)

This code estimates the vector of best-fit concentrations given the simulated concentration-normalized currents and the experimental data.


<b>obj_func_min_conc_fit_CSWV</b>(c,phi,dataInet,dataIstd)
- c = vector of best-fit concentrations to be estimated
- phi = same as in voltammogram_conc_fitting_protocol
- dataInet = same as in voltammogram_conc_fitting_protocol
- dataIstd = same as in voltammogram_conc_fitting_protocol
- return = f (the value of the objective function to be minimized)

This program generates the value of the objective function to be minimized, which is iterated over multiple times to estimate the vector of best-fit concentrations.
