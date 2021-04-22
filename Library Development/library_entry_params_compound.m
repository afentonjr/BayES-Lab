function [paramsfound,paramslb,paramsub,ETM,PETM,dataEstep,dataInet,inetsim,Eref,simtime]=library_entry_params_compound(data,paramsin,paramguess,conc,ninitguess,bsubtract,data_ref,data_bgref,prior)
%% Header
% library_entry_params_compound is a function that estimates the parameters
% and the electron transfer mechanism of a compound given an experimental
% voltammograms at known concentrations. At this current iteration, this
% function can only distinguish between two electron transfer mechanisms -
% a reversible and a quasireversible electron transfer, although additional electron
% transfer mechanisms can be added later on. The mechanisms, along with the
% estimated parameters, can be used in a Bayesian framework to determine the identity of compounds in
% solution (performed in a separate function). library_entry_params_compound feeds into the function
% "data_processing_library_find", "biccompx", and
% "voltammogram_params_fitting_protocol", which
% determines the parameters that generate the best fit, given the data,
% by maximizing the likelihood function (i.e., minimizing the exponential
% argument in the joint normal probability distribution function).
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-20
% Version: 1.0
%
%% Inputs:
%
% data -        the raw experimental data (a CSWV, with currents and potentials.
%               Potentials are in V vs. a reference redox event)
%
%               Each dataset has four columns:
%
%               1) Baseline potential (the average
%               potential of the two pulses). Since this is an average of the raw readings, the second half
%               of the column should read "NaN". Units of V vs. Ref.
%
%               2) Difference current (the difference
%               between the currents of the forward and reverse pulses at a given
%               baseline potential, with each measurement being averaged over the last
%               30% of the current pulse). Similar to the baseline potential, the second
%               half of the column should read "NaN". Should be reported in microamperes
%               (the convention for difference current used in EC-Lab).
%               Units of microA.
%
%               3) Raw potentials - the
%               potential of each forward and reverse pulse. There should be twice as
%               many entries as for baseline potential, and this column is related to
%               column 1) via the formula "Baseline_Potential(i) = mean([True_potential(2*i-1),
%               True_potential(2*i)])". Units of V vs. reference.
%
%               4) Raw current - the current,
%               averaged over the last 30% of each pulse. This current is related to
%               column 2) via the formula "Difference_Current(i) = True_Current(2*i) -
%               True_Current(2*i-1)". Should be reported in mA (the convention for
%               regular, non-difference currents in EC-Lab).
%
%               Multiple datasets should be included. The average dataset will be taken
%               to estimate the parameters.
%
%               The raw data must be pre-processed before being added to this code.
%               The peaks from different trials must be aligned. Potential bounds
%               for different trials may result in different indices for the peak
%               current, and the length of the potential and current vectors may differ
%               between trials as well. At this iteration of code, this must be done manually. If the
%               vector lengths are not equal, the code will enter debug mode. There, the
%               user needs to investigate which points to exclude from the datasets with
%               more points. The extreme potentials from these longer datasets will need
%               to be eliminated, such that the baseline potentials are aligned (e.g., the
%               baseline potentials for all trials range from -0.3 V vs. Fc to +0.7 V vs.
%               Fc. In this dataset, if one of the datasets start at -0.4 V vs. Fc, this would not be acceptable).
%               Additional pre-processing (iR correction, referencing the potential to a
%               well-known redox couple such as ferrocene, and concentration
%               normalization) are all automatically performed in this code.
%
%
% paramsin -    relevant parameters for experiment. See the SI for an example
%               CSWV waveform. paramsin must be a 1x15 cell.
%               The order is:
%               {1}  electrode radius (in mm)
%               {2}  CSWV step height (in mV)
%               {3}  CSWV pulse height (in mV)
%               {4}  CSWV pulse duration (in ms)
%               {5}  temperature (in K)
%               {6}  number of simulated points per second (units of 1/s)
%               {7}  nu, the scan rate (in mV/s)
%               {8)  R, a vector of solution resistances determined by the technique of choice (in
%                    Ohms - should be equal to the number of datasets)
%               {9}  comp, a vector of the iR compenstation automatically implemented by the program of
%                    choice (e.g., 0.85). Unitless.
%               {10} Rbg, a vector of solution resistances for the
%                    background currents (otherwise the same as R; units in Ohms)
%               {11} compbg, similar to comp but for the background
%                    currents. Unitless
%               {12} Rref, a vetor of solution resistances for the trials used for potential referencing
%                    (e.g., experiments with ferrocene added to the
%                    solution). Units of Ohms.
%               {13} compref, the compensation for the reference trials.
%                    Unitless.
%               {14} Rbgref, a vector of solution resistances for the trials used for potential referencing of the
%                    background (e.g., experiments with ferrocene added to
%                    the solution). Units of Ohms
%               {15} compbgref, the compensation for the background
%                    reference trials (unitless)
%
%
% paramsguess - a vector of guesses for parameter optimization. Must be a
%               1x5 vector. Order is:
%               (1) redox potential (V vs. ref)
%               (2) D_o (cm^2/s)
%               (3) D_r (cm^2/s)
%               (4) ko (cm/s)
%               (5) alpha. Here, "o" refers to the oxidized form of the redox
%                   couple, and "r" refers to the reduced form of the redox
%                   couple. Unitless
%
%
% conc -        known concentration (a vector), necessary for correctly determining diffusion
%               coefficients (in mol/L). The concentration for each trial should be
%               reported (e.g., if twelve CSWV voltammograms are being studied, twelve
%               concentrations must be reported).
%
%
% ninitguess -  number of randomized initial guesses (not including the
%               provided initial guess). Multiple initial guesses
%               increases the likelihood that the true global minimum is found. The code
%               is structured such that multiple trials can be run in
%               parallel. Unitless
%
%
% bsubtract -   the background current that will be subtracted from the
%               voltammogram to reduce the effects of noise. The same format and units
%               as "data". The potential range of bsubtract should be wider than that of
%               any dataset so that all datasets may be background-subtracted. In units
%               of V vs. reference (Ref for short).
%
%
% data_ref -    the data that has the compound of interest and the compound to
%               reference the potential (e.g., unsubstituted phenothiazine with
%               ferrocene). The same format and units as "data". This data is used to
%               reference the potential to that of the desired redox event. If no
%               potential adjustment is desired, input as "0" (a scalar).
%               Units: V vs. Ref
%
%
% data_bgref -  the data with only the compound to reference the potential
%               (e.g., only ferrocene). The same format and units as "data". This data is
%               used to reference the potential for the background currents. If no
%               potential adjustment is desired, input as "0" (a scalar).
%               Units: V vs. Ref
%
%
% prior -       the prior probability assigned to each electron transfer mechanism.
%               Should be a 1x2 vector (as only two electron transfer mechanisms are
%               currently considered). The sum of the vector should be 1 (i.e., normalized
%               probabilities). The first entry is the probability of a reversible
%               electron transfer, while the second entry is the probability of a
%               quasireversible electron transfer (i.e., kinetic limitations play a role
%               in the rate of the electron transfer). Is unitless.
%
%
%% Outputs:
%
% paramsfound - a cell with the parameters for each ETM. The last entry
%               denotes the estimated electron transfer mechanism. Currently, {1} is
%               a reversible electron transfer, while (2) is a quasireversible electron
%               transfer. Units vary depending on the parameter that is output. The order
%               is (1) formal potential (V vs. ref), (2) oxidized diffusion coefficient
%               (cm^2/s), (3) diffusion coefficient of reduced species (cm^2/s), (4)
%               heterogeneous rate constant (cm/s, if applicable), and (5) transfer
%               coefficient (unitless, if applicable).
%
% ETM -         the predicted electron transfer mechanism. "1" is reversible and
%               "2" is quasireversible. Unitless
%
% PETM -        the probability associated with the electron transfer. Determined
%               using the BIC, the Bayes Factor, and Bayes Rule. Unitless
%
% dataEstep -   the CSWV baseline potentials (data). Post-processed (although I
%               believe processing only affects the current). In units of V vs. ref
%
% dataInet -    the CSWV difference current (data). Post-processed. In units
%               of A.
%
% inetsim -     the CSWV difference current (simulated). In units of A.
%
% Eref -        the redox potential of the reference redox couple with respect to
%               the experimental potential (V vs. Ref). Another way to look at it
%               is the amount the potentials are adjusted by.
%
% simtime -     the simulation time. In units of seconds.


%% Start timer (for internal consumption)
tic

%% Convert from structure to doubles or cells (uncomment only if running on terminal for HPC - will need to change for different compounds)
% bsubtract=bsubtract.bsubtract;
% conc=conc.conc_inPT;
% data=data.dataPT;
% data_bgref=data_bgref.dataFc_trunc;
% data_ref=data_ref.dataPTref;
% paramguess=paramguess.params_guessPT;
% paramsin=paramsin.params_inPT;
% prior=prior.priorPT;

%% Process data beforehand
% Pre-process data, including iR correction, potential referencing,
% normalizing current signals to a single concentration, and background
% subtraction
[dataV,dataInet,dataIstd,dataEstep,Eref]=data_processing_library_find(data,bsubtract,data_ref,data_bgref,0,1,paramsin,conc);

%% Guess the params and ET mechanism
% Fit a modeled curve to experimental data for different electron transfer
% mechanisms. The estimated parameters (params), the simulated current for a given
% electron transfer mechanism (inetsim), the lower-bound estimate of the parameters
% (paramslb), and the upper-bound estimate (paramsub) are generated.
% Each cell of params is the parameters for a given electron transfer model
% each column of inetsim is a different electron transfer mechanism
% esimate. The first cell / column is for the reversible electron transfer
% mechanism, and the second cell / column is for the quasireversible electron
% transfer mechanism.
% So far, the data can only be pre-processed for potential sweeps where the
% initial sweep is from negative to positive potentials (which is the case
% on the studies with phenothiazine, where the reduced form is initially
% present). The code will be made more robust in the future for studies
% where the initial sweep is from positive to negative potentials.
[params,inetsim,paramslb,paramsub]=voltammogram_params_fitting_protocol(dataV,dataInet,dataIstd,ninitguess,paramsin,paramguess,1e-3); % Assumes a normalized concentration (1 mM)

% This function estimates the BIC for each electron transfer mechanism,
% which is used in turn to estimate the Bayes Factor
BIC=biccompx(dataInet,dataIstd,inetsim,data); % Short for BIC of component "x"
BICmin=min(BIC); % Find the minimum BIC for comparison
BF=exp((BICmin-BIC)./2); % Bayes Factor estimate
P=BF.*prior./sum(BF.*prior); % Probability of each electron transfer mechanism
[PETM,indETM]=max(P); % Find the type of electron transfer mechanism with
% the largest probability, and the respective probability
paramsfound=params; % Output the params
ETM=indETM; % Both ETM and index to use to find true params

%% Save outputs (uncomment if running on HPC)

% paramsfound,paramslb,paramsub,ETM,PETM,dataEstep,dataInet,inetsim,Eref,simtime

% save('2021-03-29-Estimated-PT-Params-100ig-1000pps.mat','paramsfound')
% save('2021-03-29-Estimated-PT-Params-LB-100ig-1000pps.mat','paramslb')
% save('2021-03-29-Estimated-PT-Params-UB-100ig-1000pps.mat','paramsub')
% save('2021-03-29-Estimated-PT-ETM-100ig-1000pps.mat','ETM')
% save('2021-03-29-Estimated-PT-ETM-Prob-100ig-1000pps.mat','PETM')
% save('2021-03-29-Estimated-PT-Potential-Data-100ig-1000pps.mat','dataEstep')
% save('2021-03-29-Estimated-PT-Current-Data-100ig-1000pps.mat','dataInet')
% save('2021-03-29-Estimated-PT-Current-Simulated-100ig-1000pps.mat','inetsim')
% save('2021-03-29-Estimated-PT-Reference-Potential-100ig-1000pps.mat','Eref')

%% Stop timer
simtime=toc;
% save('2021-03-29-Estimated-PT-Sim-Time-100ig-1000pps.mat','simtime')
end

function [dataV,dataInet,Idatastd,dataEstep,Eref]=data_processing_library_find(data,bsubtract,data_ref,data_bgref,iter,option,paramsin,conc)
%% Header
% This function pre-processes raw experimental data by first 1) iR
% correcting it, 2) referencing the potentials to that of ferrocene, and 3)
% normalizing the data so the the concentration of all the trials is the
% same. Once this adjustment is performed, then the data is averaged, the
% background current subtracted, and the standard deviation recorded for
% later in the code. All inputs come from
% "library_entry_params_compound.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
%
% data - raw experimental data. See the parent function for details
% on the formatting of this data.
%
% bsubtract - the raw background current. See the parent function for details
% on the formatting of this data.
%
% data_ref - raw experimental data, including the reference compound, for
% potential referencing (e.g., ferrocene) of the data. See the parent function for details
% on the formatting of this data.
%
% data_bgref - raw experimental data, including the reference compound, for
% potential referencing (e.g., ferrocene) of the background current.
% See the parent function for details on the formatting of this data.
%
% iter - a counter for the number of trials to be considered. If iter is
% set to 0, all trials will be considered. If only one particular trial of
% data is to be considered, then iter should be set to the desired trial
% number (ordered the same as the trials in the data - 1, 2, etc...).
%
% option - allows for background correction. If option=0, no background
% correction is performed. If option is any other value, background
% correction is performed.
%
% paramsin - the input parameters from the parent function, which contains
% more details on the formatting of these parameters.
%
% conc - the concentrations of the signals. See the parent function for
% details on the formatting of these parameters.
%
%
%% Outputs:
%
% dataV - the raw potential values (i.e., the values of the forward and
% reverse pulses), adjusted via iR compensation and referencing vs. the
% redox event of interest. In V vs. ref.
%
% dataInet - the processed difference current, normalized by the
% concentration. If iter is 0, it is the average over all trials; if iter
% is not zero, than it is that number of trial. In A.
%
% Idatastd - the standard deviation in the difference current at a given
% potential point over all trials. The effects of the variance in the
% background current are also accounted for by evaluating the covariance
% matrix of the two currents. In A.
%
% dataEstep - the baseline potential, adusted via iR compensation and
% referencing vs. the ferrocene redox event. In V vs. ref.
%
% Eref - the redox potential of the reference redox couple with respect to
%        the experimental potential (V vs. Ref). Another way to look at it
%        is the amount the potentials are adjusted by.


%% Initialization
% Set a bound so that the index does not exceed the number of array
% elements when accounting for the variance in the background data.
% This bound needs to be set in place because the number of background
% trials may not be the same as the number of data points.
minindex=min([size(data,2)/4,size(bsubtract,2)/2]);

% Extract the dataV, dataInet, and dataEstep vectors
for i=1:size(data,2)/4
    dataInet(:,i)=data(:,4.*i-2)./1e6; % Convert from mu_A to A
    dataV(:,i)=data(:,4.*i-1);
    dataEstep(:,i)=data(:,4.*i-3);
    dataI(:,i)=data(:,4*i)./1000; % Convert from mA to A
end

% Get rid of NaN or 0 - only the first half of entries for Estep
% and Inet are significant
dataInet=dataInet(1:end/2,:);
dataEstep=dataEstep(1:end/2,:);

% Extract resistances and compensations
R=paramsin{8};
comp=paramsin{9};
Rbg=paramsin{10};
compbg=paramsin{11};
Rref=paramsin{12};
compref=paramsin{13};
Rbgref=paramsin{14};
compbgref=paramsin{15};

%% IR Correction on Raw Data
for i=1:size(dataEstep,2)
    [dataEstep(:,i),dataV(:,i)]=iR_correction(dataV(:,i),dataI(:,i),R(i),comp(i));
end

%% Potential Adjustment on Raw Data
% In non-aqueous solutions, ferrocene or cobaltocene can be
% used as a reference. If cobaltocene is used, this, in turn, can be
% referenced to ferrocene so that all potentials are referenced to ferrocene.
[dataEstep,dataV,Eref]=reference_correction(dataEstep,dataInet,dataV,data_ref,Rref,compref,1);

%% Aligning Peaks for Analysis
% It is challenging to code a protocol to robustly align peaks, so we leave
% this to the user to manually perform. Guidelines are mentioned in the header of the
% parent code. If the peaks aren't aligned, the code will enter debug mode,
% encouraging the reader to explore the baseline potential waveform and make
% adjustments to the input data to align the peaks. The code will then
% error out.

TF=isnan(dataEstep); % If the input data is not the same length (i.e,. if any "NaN" remains)
if sum(sum(TF))>0
    disp('Check potential waveform and manually truncate it!')
    keyboard
    error('Adjust the input data, and then run the code again.')
else
end

%% Normalize the Currents by Concentration
for i=1:size(dataInet,2)
    dataInet(:,i)=dataInet(:,i)*1e-3./conc(i); % All currents are normalized to 1 mM
    dataI(:,i)=dataI(:,i).*1e-3./conc(i);
end

%% Background Current Processing
% Average the baseline potentials
dataEstep=mean(dataEstep,2);

% Process background subtraction vectors
for i=1:size(bsubtract,2)/4
    bgP(:,i)=bsubtract(:,4*i-3);
    bgI(:,i)=bsubtract(:,4*i-2)./1e6; % Convert from mu_A to A
    bgPraw(:,i)=bsubtract(:,4*i-1);
    bgIraw(:,i)=bsubtract(:,4*i)./1000; % Convert from mA to A
end

% Get rid of NaN or 0 - only the first half of entries for Estep
% and Inet are significant
bgP=bgP(1:end/2,:);
bgI=bgI(1:end/2,:);

%% iR-Correct Background Data
for i=1:size(bgP,2)
    [bgP(:,i),bgPraw(:,i)]=iR_correction(bgPraw(:,i),bgIraw(:,i),Rbg(i),compbg(i));
end

%% Adjust the potential to that vs. the desired reference potential
bgP=reference_correction(bgP,bgI,bgPraw,data_bgref,Rbgref,compbgref,0);

%% Average the background currents
bgPmean=mean(bgP,2);
bgImean=mean(bgI,2);
bgIstd=std(bgI,0,2);

%% Aligning the Background Current with the Raw Data
% The background data should be wider than any experimental data. These
% lines will select which sets of background data need to be extracted for
% subtraction
distmin=abs(min(dataEstep)-bgPmean(1:(length(bgPmean)+1)/2)); % This is a
% vector of distances between the minimum baseline potential of the data
% and the potential value of the background current in the forward sweep
[~,imin]=min(distmin); % Find the index that has the minimum distance
% between the baseline potential of the background data and the initial
% baseline potential of the data
imax=imin+(length(dataEstep)-1)/2; % The -1 os present because the first entry is excluded (it's already counted for in imin)

% Extract the background current to be adjusted
Iadjust=[bgImean(imin:imax);bgImean(end-(imax-2):end-(imin-1))];
Ibgexp=[bgI(imin:imax,:);bgI(end-(imax-2):end-(imin-1),:)];
% imax-2, rather than imax-1, is used because the
% Biologic potentiostat reaches the max baseline potential on the forward
% sweep. It does not repeat the max potential on the reverse sweep. For
% example, the baseline potential waveform is of the form 0.5 V, 0.6 V,
% 0.7 V, 0.6 V, etc...; instead of 0.5 V, 0.6 V, 0.7 V, 0.7 V, 0.6 V,...
% However, if the max potential is repeated, then the imax-2 can be
% replaced with imax-1, as done below.
if length(Iadjust)==length(dataEstep)-1
    Iadjust=[bgImean(imin:imax);bgImean(end-(imax-1):end-(imin-1))];
    Ibgexp=[bgI(imin:imax,:);bgI(end-(imax-1):end-(imin-1),:)];
end

%% Standard Deviation Calculation
% Take the standard deviation of the data
dataIstd=std(dataInet,0,2);

if option==0 % Option to not BG correct
    Idatastd=dataIstd;
else % BG Correct
    A=[1,-1]; % The raw current is being added (+1), and the background current is being subtracted (-1)
    for i=1:size(Ibgexp,1)
        V=cov(dataInet(i,1:minindex),Ibgexp(i,1:minindex)); % Create a covariance matrix between the raw current and background current
        Idatastd(i)=A*V*A'; % Calculate the covariance
        if Idatastd(i)<0
            error('Unphysical system - variance is negative')
        end
        Idatastd(i)=sqrt(Idatastd(i)); % Take the square root to get the std dev
    end
    Idatastd=Idatastd'; % Convert to a column vector
    dataInet=dataInet-Iadjust;
end

%% Final Processing
% Calculate mean absolute potential, difference current, and baseline
% potential
if iter==0
    dataInet=mean(dataInet,2); % Average the trials
elseif iter>size(dataInet,2)
    error('Iter is larger than the number of trials. Make iter smaller.')
else
    dataInet=dataInet(:,iter); % Test an individual trial
end
dataV=mean(dataV,2); % Average dataV
%%
end

function [dataEstep,dataVnew]=iR_correction(dataV,dataI,R,comp)
%% Header
% This function iR-corrects both the raw potential and the baseline
% potential, given the true current, the resistance, and the compensation.
% All inputs are directly fed from "data_processing_library_find.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
% dataV - the raw (not baseline) potential values, in V vs. ref.
%
% dataI - the raw (not difference) currents, in A.
%
% R - the solution resistance, in Ohms.
%
% comp - the compensation (in a decimal form between 0 and 1). 1-comp is the uncompensated (e.g., at
% 85% ZIR, comp=0.85, 1-comp=0.15)
%
%
%% Outputs:
% dataEstep - the iR-corrected baseline potential, in V vs. ref
%
% dataVnew - the iR-corrected raw potential, in V vs. ref
%
dataVnew=dataV-dataI*(1-comp)*R; % Perform a ZIR correction on the raw potential
for i=1:length(dataV)/2
    dataEstep(i)=(dataV(2*i-1)-dataI(2*i-1)*(1-comp)*R+dataV(2*i)-dataI(2*i)*(1-comp)*R)/2; % Compensating the baseline potential
end
dataEstep=dataEstep';
end

function [dataEstep,dataVnew,Eref]=reference_correction(dataEstep_comp,dataInet_comp,dataV,data_ref,Rref,comp,bg)
%% Header
% This function takes in potentials measured against some pseudoreference
% (e.g., Ag/Ag+) and adjusts the potential to some desired reference redox
% event (e.g., ferrocene) using data that has the redox event of interest
% (if any) and the reference redox event. This code can both change the
% reference potentials for the experimental data and the background data.
% All inputs are fed from "data_processing_library_find.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
% dataEstep_comp - the baseline potential for the compound (if any; hence "comp")
% of interest. No reference compound is present. In units of V vs. ref.
%
% dataInet_comp - the difference current for the compound (if any) of interest.
% No reference compound is present. In units of V vs. ref.
%
% dataV - the raw potential including the compound (if any) of interest.
% No reference compound is present. In units of V vs. ref.
%
% data_ref - the full matrix (i.e., currents and potentials) of all data containing the reference compound.
% This full matrix will be further processed in the code below.
% Directly fed from "data_processing_library_find.m", which in turn is
% fed from "library_entry_params_compound.m", where further details may be
% found.
%
% Rref - the resistance of the solution. In units of Ohms.
%
% comp - the compensation (in a decimal form between 0 and 1). 1-comp is the uncompensated (e.g., at
% 85% ZIR, comp=0.85, 1-comp=0.15)
%
% bg - an index determining whether the data whose potential is being
% corrected is the background current. Input 0 if the data is the background
% data, and input another value if the data is the experimental data (i.e.,
% containing the compound of interest).
%
%
%% Outputs:
%
% dataEstep - the new baseline potential, referenced to the desired redox
% event of interest. In V vs. ref
%
% dataVnew - the new raw potential, referenced to the desired redox event
% of interest. In V vs. ref
%
% Eref - the redox potential of the reference redox couple with respect to
%        the experimental potential (V vs. Ref). Another way to look at it
%        is the amount the potentials are adjusted by.

%% Main Body
% First, find the maximum potentials of the peak from the compound being
% tested
% There should only be two peaks of significance in the data
% Loop over all datasets here
for i=1:size(dataEstep_comp,2)
    [~,ind_comp]=findpeaks(abs(dataInet_comp(:,i)),'SortStr','descend','NPeaks',2);
    E0_comp(i)=mean([dataEstep_comp(ind_comp(1),i),dataEstep_comp(ind_comp(2),i)]);
end

%% Extract the potentials and currents from the reference data if it exists
if max(max(data_ref))==0 % If no background-correction data is put in, make no adjustment
    warning('No adjustment will be made to another reference potential.')
    Eref=0;
    dataEstep=dataEstep_comp;
    dataVnew=dataV;
else
    for i=1:size(data_ref,2)/4
        refP(:,i)=data_ref(:,4*i-3);
        refI(:,i)=data_ref(:,4*i-2)./1e6; % Convert from mu_A to A
        refPraw(:,i)=data_ref(:,4*i-1);
        refIraw(:,i)=data_ref(:,4*i)./1000; % Convert from mA to A
    end
    
    % Get rid of NaN or 0 - only the first half of entries for Estep
    % and Inet are significant
    refP=refP(1:end/2,:);
    refI=refI(1:end/2,:);
    
    %% iR-Correct the reference data
    for i=1:size(refP,2)
        [refP(:,i),refPraw(:,i)]=iR_correction(refPraw(:,i),refIraw(:,i),Rref(i),comp(i));
    end
    
    %% Integrate the reference data
    % Loop over all reference datasets here
    for j=1:size(refP,2)
        [~,ind_ref,~,p_ref]=findpeaks(abs(refI(:,j)));
        if bg==0 % Background current potential adjustment
            [~,ind_ref,~,~]=findpeaks(abs(refI(:,j)),'SortStr','descend','NPeaks',2); % Limit to two peaks if the background current is being corrected.
            Eref=mean([refP(ind_ref(1),j),refP(ind_ref(2),j)]);
            dataEstep(:,j)=dataEstep_comp(:,j)-Eref; % Adjust potentials
            dataVnew(:,j)=dataV(:,j)-Eref;
        else % Non-background current potential adjustment
            % Eliminate peaks with a prominence less than 10 times that of the peak
            % with the largest prominence. This is an emperical threshold
            % that can be adjusted if desired.
            elim_vec=0;
            for i=1:length(p_ref)
                if p_ref(i)<0.1*max(p_ref)
                    elim_vec(i)=1;
                else
                    elim_vec(i)=0;
                end
            end
            ind_ref(elim_vec==1)=[];
            p_ref(elim_vec==1)=[];
            
            if length(p_ref)==2 % This would only occur if the reference
                % redox event of interest has a convoluted peak with that
                % from the compound whose properties are being extracted
                warning('No adjustment will be made to another reference potential.')
                Eref=0;
                dataEstep=dataEstep_comp;
                dataVnew=dataV;
            elseif length(p_ref)==4
                ind_ref=sort(ind_ref); % Sort the indices in ascending order - one extrema (either maxima or minima) will be the first two, the other will be the second two
                Ediffj=0.5.*(abs(refP(ind_ref(1),j)-refP(ind_ref(2),j))+...
                    abs(refP(ind_ref(3),j)-refP(ind_ref(4),j))); % Find the difference in the redox potentials
                Esign=mean([refP(ind_ref(1),j),refP(ind_ref(2),j),refP(ind_ref(3),j),refP(ind_ref(4),j)]); % Find the average of all peaks
                if Esign-E0_comp>0 % If the average of the two peaks is greater than the average of the peak from the compound of interest
                    Ediff=-Ediffj; % Switch the sign of Ediff
                else
                    Ediff=Ediffj;
                end
                Eref=E0_comp(j)-Ediff; % Calculate the reference potential
                dataEstep(:,j)=dataEstep_comp(:,j)-Eref; % Adjust potentials
                dataVnew(:,j)=dataV(:,j)-Eref;
            else
                error('Too many peaks identified')
            end
        end
    end
end
end

function [params,inetsim,lowerbound,upperbound]=voltammogram_params_fitting_protocol(dataV,dataInet,dataIstd,ninitguess,paramsin,paramguess,conc)
%% Header
% This function takes the input experimental data, already pre-processed,
% and estimates the parameters, the optimal simulated difference current,
% the lower bound estimate of the parameters (determined via an emperical
% threshold), and an upper bound estimate, also set to the same emperical
% threshold. All inputs are directly fed from
% "library_entry_params_compound.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
%
% dataV - the raw, processed potential data, in units of V vs. ref.
%
% dataInet - the difference current that has been processed. In units of A.
%
% dataIstd - the standard deviation of the processed difference current,
% experimentally determined. This also accounts for variations in the
% background current. In units of A.
%
% ninitguess - the number of randomized initial guesses (not including the
% initial guess provided by the user).
%
% paramsin - the input parameters.
%
% paramsguess - the guess for the specific parameters to be estimated.
%
% conc - the concentration of the system. This should be scalar, as the
% difference currents from all the trials should be normalized to a
% singular concentration.
%
%
%% Outputs:
%
% params - the estimated parameters. Is a cell. The first cell
% is for a reversible electron transfer. The second cell is a
% quasireversible electron transfer. Units vary depending on the parameter.
%
% inetsim - the simulated current that best matches the experimental data
% for a given electron transfer mechanism. The first column is for a
% reversible electron transfer mechanism. The second column is for a
% quasireversible electron transfer mechanism. In units of A.
%
% lowerbound - a lower bound of the estimates for each of the parameters.
% This is determined by an emperical standard. Currently, for a given
% parameter of interest, when all other parameters are kept at their
% optimal value, the value of the lower bound is that which results in a
% liklihood half that of the maximum likelihood. Same units as in "params".
%
% upperbound - an upper bound estimate for each of the parameters, using
% the same criteria as that for the lower bound estimates.

%% Import parameters
re=paramsin{1};
dE=paramsin{2}/1000; % Convert dE, Esw, and tpulse to SI units here (from mV, mV, ms to V, V, s, respectively)
Esw=paramsin{3}/1000;
tpulse=paramsin{4}/1000;
T=paramsin{5};
ptspersec=paramsin{6};
nu=paramsin{7};
dnu=dE/(ptspersec*2*tpulse);
%% Create lower and upper bound estimates for fmincon
Efguess=paramguess(1); % In V vs. ref
Doguess=paramguess(2); % cm^2/s
Drguess=paramguess(3); % cm^2/s
Depslow=0.5.*min(Doguess);
Depshigh=max(Doguess);
% The order of parameters for the bounds: formal potential (Ef, V vs. ref), diffusion
% coefficient of the oxidized species (Do, cm/s), diffusion coefficient of the reduced
% species (Dr, cm/s), heterogeneous rate constant (ko, cm/s), and the electron transfer
% coefficient (alpha, unitless)
lowerb=[(min(Efguess)-0.03),(min(Doguess)-Depslow),(min(Drguess)-Depslow),1e-5,0]; % Lower bound for the allowable range used by fmincon
upperb=[(max(Efguess)+0.03),(max(Doguess)+Depshigh),(max(Drguess)+Depshigh),1e1,1]; % Upper bound
for i=1:ninitguess+1 % The more initial guesses are employed, the better the fit will be, given that this is a 3-5 dimensional optimization
    % Parameter guess rng
    pguess(i,:)=lowerb+rand(size(lowerb)).*(upperb-lowerb);
    % Since ko is better suited for a log scale, the rng is slightly
    % modified for this value
    randko=ceil(6*rand); % Determine the order of magnitude as a random integer
    if randko==1 % If the number is 1
        pguess(i,4)=1e-5+rand.*(1e-4-1e-5); % ko guess is in the range of 1e-5 to 1e-4 cm/s
    elseif randko==2
        pguess(i,4)=1e-4+rand.*(1e-3-1e-4);
    elseif randko==3
        pguess(i,4)=1e-3+rand.*(1e-2-1e-3);
    elseif randko==4
        pguess(i,4)=1e-2+rand.*(1e-1-1e-2);
    elseif randko==5
        pguess(i,4)=1e-1+rand.*(1e0-1e-1);
    elseif randko==6
        pguess(i,4)=1e0+rand.*(1e1-1e0);
    end
end
for j=1:2 % Test over all possible ETM models
    if j==1 % Reversible 1 e- model
        pguess0=pguess(:,1:3); % Only the first three parameters (Eo, Do, Dr) matter
        lowerb0=lowerb(1:3);
        upperb0=upperb(1:3);
        cc=zeros(1,length(upperb0)); % Preallocate vectors for the parameters
        dd=zeros(1,length(upperb0));
        cguessend=paramguess(1:3); % The final initial guess vector is the user-defined initial guess (is a row vector)
    elseif j==2 % Quasireversible 1 e- model
        pguess0=pguess; % All five parameters matter here
        lowerb0=lowerb;
        upperb0=upperb;
        cc=zeros(1,length(upperb0));
        dd=zeros(1,length(upperb0));
        cguessend=paramguess;
    end
    parfor i=1:ninitguess+1 % Parallelizing
        cguess=pguess0(i,:); % The randomized initial guess
        if i==ninitguess+1
            cguess=cguessend; % The user-input initial guess
        end
        params_fit=[conc,re,dE,Esw,tpulse,T,dnu,nu];
        % Different algorithms are used for a given initial guess to better
        % fit the parameters to data. The fit with the maximum likelihood
        % is accepted
        options=optimoptions(@fmincon,'Algorithm','interior-point','FiniteDifferenceType','central','Display','Iter','StepTolerance',1e-12,'TolFun',1e-12,...
            'MaxFunctionEvaluations',500,'SubproblemAlgorithm','cg');
        [cc(i,:),fval1(i)]=fmincon(@(c)obj_func_min_param_fit_CSWV(c,params_fit,dataInet,j,dataV,dataIstd),...
            cguess,[],[],[],[],lowerb0,upperb0,[],options);
        options1=optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','central','Display','Iter','StepTolerance',1e-12,'TolFun',1e-12,...
            'MaxFunctionEvaluations',500);
        [dd(i,:),fval2(i)]=fmincon(@(c)obj_func_min_param_fit_CSWV(c,params_fit,dataInet,j,dataV,dataIstd),...
            cguess,[],[],[],[],lowerb0,upperb0,[],options1);
        subcycleperc=i/(ninitguess+1) % Display a counter to track progress
    end
    % Evaluate the minimum of all the fits for all initial guesses
    [min1,ind1]=min(fval1);
    [min2,ind2]=min(fval2);
    if min1>min2 % Compare the fits for a given method
        params{1,j}=dd(ind2,:); % Accept the grand minimum only
        options2=optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','central','Display','Iter','StepTolerance',1e-12,'TolFun',1e-12,...
            'MaxFunctionEvaluations',500); % For the determination of the bounds, use the set of options that found the best fit (to save time)
    elseif min1<=min2
        params{1,j}=cc(ind1,:);
        options2=optimoptions(@fmincon,'Algorithm','interior-point','FiniteDifferenceType','central','Display','Iter','StepTolerance',1e-12,'TolFun',1e-12,...
            'MaxFunctionEvaluations',500,'SubproblemAlgorithm','cg');
    end
    params_fit=[conc,re,dE,Esw,tpulse,T,dnu,nu]; % For importation into the subsequent for loop
    parfor q=1:length(params{1,j}) % Parallelize over the lower bounds first
        lb=lowerb(q); % The lower bound is the same as for the parameter estimation protocol
        ub=params{1,j}(q); % The upper bound is the estimate of the parameters from the previous fitting protocol
        lbfind(q)=fmincon(@(c)obj_func_find_param_bounds(c,params_fit,dataInet,j,dataV,dataIstd,min([min1,min2]),params{1,j},q),...
            params{1,j}(q),[],[],[],[],lb,ub,[],options2);
    end
    lowerbound{1,j}=lbfind;
    parfor q=1:length(params{1,j}) % Then parallelize over the upper bounds
        lb=params{1,j}(q); % The lower bound here is the estimate of the parameters from the previous fitting protocol
        ub=upperb(q); % Keep the same upper bound
        ubfind(q)=fmincon(@(c)obj_func_find_param_bounds(c,params_fit,dataInet,j,dataV,dataIstd,min([min1,min2]),params{1,j},q),...
            params{1,j}(q),[],[],[],[],lb,ub,[],options2);
    end
    upperbound{1,j}=ubfind;
    if j==1
        Ef=params{1,j}(1); % Assign output values to meaningful parameters
        Do=params{1,j}(2);
        Dr=params{1,j}(3);
        ko=1;
        alpha=1;
    elseif j==2
        Ef=params{1,j}(1);
        Do=params{1,j}(2);
        Dr=params{1,j}(3);
        ko=params{1,j}(4);
        alpha=params{1,j}(5);
    end
    params_sim=[conc,Dr,Do,re,dE,Esw,tpulse,T,Ef,alpha,ko,dnu,nu]; % Assign the parameter vector to simulate the best-fit voltammogram, for output
    if dataV(3)>dataV(1) % If the initial sweep is oxidative
        [~,~,inetsim(:,j),~,~]=Voltammetry_Simulator_MacroDisk_2019_11_07(j,1,2,dataV,params_sim);
    elseif dataV(3)<dataV(1) % If the initial sweep is reductive
        [~,~,inetsim(:,j),~,~]=Voltammetry_Simulator_MacroDisk_2019_11_07(j,2,2,dataV,params_sim);
    end
end
end

function BIC=biccompx(dataInet,dataIstd,inetsim,data)
%% Header
% This function calculates the BIC of an electron transfer mechanism given
% the experimental data, the standard deviation, the simulated current, and
% the original input data. All inputs are fed from
% "library_entry_params_compound.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
%
% dataInet - the processed difference current, in units of A.
%
% dataIstd - the experimental standard deviation, accounting for deviations
% in the background-subtracted current, in units of A.
%
% inetsim - the simulated current that best fits the experimental data, for
% both candidate electron transfer mechanisms. In units of A.
%
% data - the original dataset. Used to estimate the number of data points used in the calculation of the BIC.
%
%
%% Output:
%
% BIC - the output BIC value for the data being analyzed
for i=1:size(inetsim,2)
    nparams=0;
    % Determine the number of parameters
    if i==1
        nparams=nparams+3; % # of params for a reversible electron transfer mechanism
    elseif i==2
        nparams=nparams+5; % # of params for a quasireversible electron transfer mechanism
    end
    Lfull=0; % Initialize the sum
    for j=1:length(dataInet)
        Lfull=Lfull+log(1/((dataIstd(j))*sqrt(2*pi)))+(-(1/(2*dataIstd(j)^2))*(dataInet(j)-inetsim(j,i))^2); % Sum the logs
    end
    BIC(i)=log((size(data,2)./4).*length(inetsim(:,i))).*nparams-2*Lfull; % The formula for BIC is k*ln(n)-2*ln(likelihood).
    % "k" is the number of parameters
    % "n" is the number of data points. There are "m" repeats of "r" points
    % each. As such, n=m*r.
end
end

function f=obj_func_min_param_fit_CSWV(c,params_fit,dataInet,model,dataV,dataIstd)
%% Header
% Objective function to find the best estimate of the parameters, given the
% experimental data. All data imported from
% "voltammogram_params_fitting_protocol.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
%
% c - vector of parameters to be estimated. Units vary depending on the
% parameter.
%
% params_fit - fixed parameters used to fit the simulated model.
%
% dataInet - the processed difference current. Units of A.
%
% model - the electron transfer model being considered.
%
% dataV - the raw processed potential data. Units of V vs. ref.
%
% dataIstd - the standard deviation of the experimental current, including
% the standard deviation from the background current. Units of A.
%
%
%%  Outputs:
%
% f - the objective function to be minimized (in this case, the argument in
% the exponent of the maximum likelihood function).
%
%
if model==1 % Reversible ET
    % Order of parameters: conc,Dr,Do,re,dE,Esw,tpulse,T,Ef,alpha,ko,dnu,nu
    params_fit=[params_fit(1),c(3),c(2),params_fit(2),...
        params_fit(3),params_fit(4),params_fit(5),params_fit(6),c(1),1,1,...
        params_fit(7),params_fit(8)];
elseif model==2
    params_fit=[params_fit(1),c(3),c(2),params_fit(2),...
        params_fit(3),params_fit(4),params_fit(5),params_fit(6),c(1),c(5),c(4),...
        params_fit(7),params_fit(8)];
end
if dataV(3)>dataV(1) % If an oxidative sweep at first
    [~,~,inet,~,~]=Voltammetry_Simulator_MacroDisk_2019_11_07(model,1,2,dataV,params_fit);
elseif dataV(3)<dataV(1) % If a reductive sweep at first
    [~,~,inet,~,~]=Voltammetry_Simulator_MacroDisk_2019_11_07(model,2,2,dataV,params_fit);
end
inet=inet'; % Make into a column vector
ivec=(dataInet-inet).^2./dataIstd.^2; % Maximize the likelihood by weighting errors with standard deviations
f=sum(ivec); % Maximizing the likelihood, which is the product of normal distributions, is the same as adding the arguments in the exponents
if f==Inf || isnan(f)==1 || isreal(f)==0 % Used to circumnaviate errors when the objective function
                                         % is initially one of these three criteria
    f=1e20; % Set f to an arbitrary very large real value
end
end

function f=obj_func_find_param_bounds(c,params_fit,dataInet,model,dataV,dataIstd,fval,params,q)
%% Header
% Objective function to find the best estimate of the lower and upper parameter bounds,
% given the experimental data and how well the best estimate of the parameters
% fit the experimental data. All data imported from
% "voltammogram_params_fitting_protocol.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2021-03-10
% Version: 1.0
%
%% Inputs:
%
% c - the lower or upper bound to be estimated. Units vary depending on the
% parameter.
%
% params_fit - fixed parameters used to fit the simulated model.
%
% dataInet - the processed difference current. Units of A.
%
% model - the electron transfer model being considered.
%
% dataV - the raw processed potential data. Units of V vs. ref.
%
% dataIstd - the standard deviation of the experimental current, including
% the standard deviation from the background current. Units of A.
%
% fval - the value of the minimized argument in the exponent. Unitless
%
% params - the best estimate of the parameters. Units vary depending on the
% parameter.
%
% q - the index which is used to determine which parameter the lower or
% upper bound is being estimated for.
%
%
%% Outputs:
%
% f - the objective function to be minimized (in this case, the difference
% of the likelihood given an estimate of the bound and half of the maximum
% likelihood).
%
if model==1 % Reversible ET
    % Order of parameters: conc,Dr,Do,re,dE,Esw,tpulse,T,Ef,alpha,ko,dnu,nu
    if q==1 % For each different parameter tested, the params_fit vector needs to change accordingly. This is for Ef
        params_fit=[params_fit(1),params(3),params(2),params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),c,1,1,...
            params_fit(7),params_fit(8)];
    elseif q==2 % This is for Do
        params_fit=[params_fit(1),params(3),c,params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),params(1),1,1,...
            params_fit(7),params_fit(8)];
    elseif q==3 % This is for Dr
        params_fit=[params_fit(1),c,params(2),params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),params(1),1,1,...
            params_fit(7),params_fit(8)];
    end
elseif model==2
    if q==1 % Ef
        params_fit=[params_fit(1),params(3),params(2),params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),c,params(5),params(4),...
            params_fit(7),params_fit(8)];
    elseif q==2 % Do
        params_fit=[params_fit(1),params(3),c,params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),params(1),params(5),params(4),...
            params_fit(7),params_fit(8)];
    elseif q==3 % Dr
        params_fit=[params_fit(1),c,params(2),params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),params(1),params(5),params(4),...
            params_fit(7),params_fit(8)];
    elseif q==4 % ko
        params_fit=[params_fit(1),params(3),params(2),params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),params(1),params(5),c,...
            params_fit(7),params_fit(8)];
    elseif q==5 % alpha
        params_fit=[params_fit(1),params(3),params(2),params_fit(2),...
            params_fit(3),params_fit(4),params_fit(5),params_fit(6),params(1),c,params(4),...
            params_fit(7),params_fit(8)];
    end
end
if dataV(3)>dataV(1) % If an oxidative sweep at first
    [~,~,inet,~,~]=Voltammetry_Simulator_MacroDisk_2019_11_07(model,1,2,dataV,params_fit);
elseif dataV(3)<dataV(1) % If a reductive sweep at first
    [~,~,inet,~,~]=Voltammetry_Simulator_MacroDisk_2019_11_07(model,2,2,dataV,params_fit);
end
inet=inet'; % Make into a column vector
ivec=(dataInet-inet).^2./dataIstd.^2;
isum=sum(ivec);
f=(isum-(fval-log(0.5))).^2; % The 0.5 is indicative that the likelihood function is halved.
% The square is minimized, which is equivalent to zeroing the objective function
if f==Inf || isnan(f)==1 || isreal(f)==0 % Used to circumnaviate errors when the objective function
                                         % is initially one of these three criteria
    f=1e20; % Set f to an arbitrary very large real value
end
end