function [P,conc0,concfinal,conclbfinal,concubfinal,Estep,inet0,inet,dataInet,Eref,simtime]=conc_compound_find(data,params,ninitguess,lib,bsubtract,data_ref,data_bgref,prior,type)
%% Header
% conc_compound_find is a function that determines the concentration of a
% number of species in solution, given a set of parameters (e.g., diffusion
% coefficients) and an electron transfer mechanism organized in a library. This is then used to
% estimate the probability of a particular compound being present in the electrolyte.
% Data pre-processing is also automated, with the exception of making
% the number of potential-current values for each trial the same
% (to ensure the same vector length).
%
% This protocol is capable of identifying compounds from only one
% voltammetric experiment. This holds an advantage over having to assess
% multiple voltammetric experiments because the electrolyte composition may
% be transient over the time scale necessary to conduct repeats.
%
% However, the voltammograms taken for background correction and that of a reference
% compound must be taken in repeats. This is because the standard deviation
% for the single trial is estimated from that of these voltammograms via
% linearization of the standard deviation as a function of current value.
% Futher details are found on the list of inputs.
% This parent function feeds into the functions "data_processing_conc_find_CV.m" or
% "data_processing_conc_find_CSWV.m", "simcurrent.m",
% "voltammogram_conc_fitting_protocol.m", and "biccompx.m".
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-20
% Version: 1.0
%
%% Inputs:
% data -       the raw experimental data (either a conventional CV or a Cyclic
%              Square Wave - CSW - Voltammogram). Each CV trial must have two
%              columns. First column: potential (in V vs. Ref). Second column:
%              current (in mA - the standard for EC-Lab). Each CSWV trial
%              (i.e., a single CSWV) must have four columns. First column: baseline
%              potential (V vs. Ref). Second column: difference current (in microA,
%              the convention for EC-Lab). Third column: true potential,
%              (V vs. Ref); averaging each set of two entries results in the baseline 
%              potential. Fourth column: true current (averaged over the last 30% of 
%              each pulse), in mA. The difference of the forward and reverse current
%              yields column 2, the difference current. NOTE: using this format, there
%              should be twice as many entries for the true potential and current as
%              there are for the baseline potential and difference current. As such,
%              there are unfilled entries in this data matrix. These should be
%              with NaN for this code to properly process them - this can be
%              done by importing the data into MATLAB using "xlsread".
%
% params -     relevant parameters for experiment, listed below as a cell.
%              Note this must be a cell, rather than a vector. The order is:
%              {1}  Vector of initial concentration guesses (should be a
%                   1 by nlib vector, where nlib is the number of entries in
%                   the library being referenced). Units of mol/L
%              {2}  The radius of the working electrode. Units of mm
%              {3}  The step height of the CSW voltammogram. Units of mV
%              {4}  The pulse height of the CSW voltammogram. Units of mV
%              {5}  The pulse duration of the CSW voltammogram. Units of ms
%              {6}  The temperature of the electrolyte. Units of K
%              {7}  The time resolution of the simulation (for CSW
%                   voltammograms only - the input time resolution is used
%                   for conventional CVs). Units of points per second
%                   (1/s).
%              {8}  The true (for a conventional CV) or effective (for a CSW
%                   voltammogram) scan rate. Units of mV/s
%              {9}  The resistance of the electrolyte in which compounds are
%                   to be identified. Should be a scalar, as only one
%                   voltammogram should be evaluated. Units of Ohms
%              {10} The compensation automatically applied by the
%                   potentiostat for the electrolyte of interest.
%                   For example, if the potentiostat is operating
%                   at 85% compensation, input the value 0.85. Unitless (a
%                   fraction)
%              {11} The resistance of the background trials being
%                   considered (for background subtraction). Should be a
%                   1 by n_bg vector, where n_bg is the number of
%                   background trials considered. Units of Ohms
%              {12} The compensation automatically applied by the
%                   potentiostat for the background trials. Should be a 1
%                   by n_bg vector, with n_bg holding the same meaning.
%                   Unitless (a fraction)
%              {13} The resistance of the electrolyte used to reference the
%                   redox potential to a reference couple (e.g., ferrocene).
%                   Should be a scalar, as only one
%                   voltammogram needs to be evaluated. Units of Ohms
%              {14} The compensation automatically applied by the
%                   potentiostat for the electrolyte used to reference the
%                   redox potential to a reference couple (e.g., ferrocene).
%                   Unitless (a fraction)
%              {15} The resistance of the trials being considered to
%                   shift the redox potential of the background subtraction
%                   to a desired redox event. Should be a
%                   1 by n_bgref vector, where n_bgref is the number of
%                   trials considered for referencing the background potentials. Units of Ohms
%              {16} The compensation automatically applied by the
%                   potentiostat for the trials used to shift the potential
%                   of the background trial. Should be a 1 by n_bgref vector,
%                   with n_bgref holding the same meaning. Unitless (a fraction)
%              {17} Concentrations of the reference active species used to
%                   adjust the potential of the background currents. Should
%                   be a 1 by n_bgref vector, with n_bgref holding the same
%                   meaning. Units of mol/L
%
% ninitguess - number of initial concentration guesses (unitless). This optimization
%              is a linear optimization with constraints, so the true minimum may be
%              rapidly found. As such, not many initial guesses are likely needed
%
% lib -        library of compounds, with relevant parameters listed. The order
%              and units are:
%              1) Formal potential Ef (V vs. Ref)
%              2) Diffusion coefficient of the oxidized species Do (cm^2/s)
%              3) Diffusion coefficient of the reduces species Dr (cm^2/s)
%              4) Heterogeneous rate constant ko (cm/s, if applicable)
%              5) Transfer coefficient in the Butler Volmer Equation (unitless, if applicable)
%
%              The last entry in each cell is the electron transfer mechanism (is unitless). If the 
%              entry value is 1, then the ET model considered is reversible one-
%              electron. If the value is 2, then the ET model is quasireversible one-
%              electron. If a reversible ET model is assumed, then ko and alpha aren't
%              listed, and only four entries are listed for the compound of interest.
%              (in this case, the code interally prevents its crashing by assigning an arbitrary
%              value of 1 to both ko and alpha).
%
% bsubtract -  a dataset of blank electrolyte (only supporting salt) used to
%              estimate the effects of the background current and minimize
%              its effects from the dataset being examined. Must have more than one trial
%              (e.g., multiple repeats), as this is used to estimate the standard
%              deviation of the single experimental trial being examined.
%              The formatting and units for each trial are the same as for
%              "data" for both a conventional CV and a CSW voltammogram.
%
% data_ref -   a single trial of the solution of interest plus that of the
%              reference redox compound. This is used to reference the potential
%              of the experimental data to that of a desired redox potential.
%              The formatting and units are the same as in "data". If data_ref is set to
%              0, no potential shift will occur
%
% data_bgref - multiple trials of an electrolyte containing only electroactive species
%              being that of the reference redox compound. This is used to reference the potential
%              of the background data to that of a desired redox potential,
%              and is also used to estimate the standard deviation of the
%              experimental current by assuming that the magnitude of the
%              experimental current and its standard deviation are linearly
%              correlated (regardless of the voltammetric experiment). Even
%              if no background correction or potential referencing is
%              taking place, this data is necessary for the above reason.
%              The formatting and units are the same as in "data".
%
% prior -      the prior probability of existence assigned to each compound. Should be a 1 x n 
%              vector, where "n" is the number of compounds in the library. Is unitless.
%
% type -       determines whether CV data or CSWV data is being analyzed (is unitless).
%              If type = 1, CV data is being analyzed.
%              If type is any other value, CSWV data is being analyzed.
%
%
%% Outputs:
% P -           probability that each compound in the library is present in the
%               electrolyte being examined. The ordering is (i) PT, (ii) MPT, (iii) EPT, (iv) iPrPT,
%               and (v) PhPT The ordering depends on the ordering of the
%               library, which can be adjusted.
%
% conc0 -       concentration vector of all components considered, in units of
%               mM (or mol/m^3). The ordering is (i) PT, (ii) MPT, (iii) EPT, (iv) iPrPT,
%               and (v) PhPT. Note that even compounds deemed not in solution
%               (i.e., assigned probabilities of less than 50 %) are listed here.
%
% concfinal -   concentration vector of all components with a probability of
%               greater than 50 % of being in solution (the percentage can be adjusted), in 
%               units of mM (or mol/m^3). The ordering is PT, MPT, EPT, iPrPT,
%               and PhPT, if it is included in the solution. For example,
%               if PT, iPrPT, and PhPT are deemed probable of being in solution, then the
%               vector would be [c_PT, c_iPrPT, c_PhPT].
%
% conclbfinal - vector of lower bound concentrations for all components with a probability of
%               greater than 50 % of being in solution (the percentage can be adjusted), in 
%               units of mM (or mol/m^3). The ordering is PT, MPT, EPT,
%               iPrPT, and PhPT, if it is included in the solution.
%               The values are determined by finding the concentration of the
%               compound of interest, less than its optimal value, that
%               halves the likelihood function while all other concentrations
%               are kept at their optimal value. This metric can naturally be
%               adjusted as desired.
%
% concubfinal - vector of upper bound concentrations for all components with a probability of
%               greater than 50 % of being in solution (the percentage can be adjusted), in 
%               units of mM (or mol/m^3). The ordering is PT, MPT, EPT,
%               iPrPT, and PhPT, , if it is included in the solution.
%               The values are determined by finding the concentration of the
%               compound of interest, greater than its optimal value, that
%               halves the likelihood function while all other concentrations
%               are kept at their optimal value. This metric can naturally be
%               adjusted as desired.
%
% Estep -       the potential (baseline potential if a CSW voltammogram is being analyzed)
%               of the best fit curve to the data, in
%               units of V vs. Ref. This potential is the potential from the data
%               once it is post-processed.
%
% inet0 -       a vector of simulated difference currents found by multiplying 
%               the concentration vector conc0 with the simulated concentration-
%               normalized currents from the library. In units of A.
%
% inet -        a vector of simulated difference currents found by multiplying 
%               the concentration vector concfinal with the simulated concentration-
%               normalized currents from the library. As such, only compounds with
%               assigned probabilities of greater than 50 % are considered in
%               this fit. In units of A.
%
% dataInet -    the post-processed current data to be fitted. This is current
%               if a conventional CV is studied, and difference current if a CSW
%               voltammogram is considered. Units are in A.
%
% Eref -        the redox potential of the reference redox couple with respect to
%               the experimental potential (V vs. Ref). Another way to look at it
%               is the amount the potentials are adjusted by.
%
% simtime -     the time for the entire simulation (s)

%% Start timer
tic % For internal consumption
%% Pre-process Data
nlib=size(lib,2); % Find the size of the library
% Data pre-processing differs depending on whether a conventional CV or a
% CSW voltammogram is being studied.
if type==1 % If a conventional CV is being studied
    [dataV,dataInet,dataIstd,bgIstd,Eref]=data_processing_conc_find_CV(data,bsubtract,data_ref,data_bgref,params,nlib);
else % If a CSW voltammogram is being studied
    [dataV,dataInet,dataIstd,dataEstep,bgIstd,Eref]=data_processing_conc_find_CSWV(data,bsubtract,data_ref,data_bgref,params,nlib); % Converts current from mA to A. 1 as the second input brings forth background correction
end
% Note that reference potential correction may be done by running
% "conc_compound_find.m" up to this point, using Eref to manually adjust
% the potential, and then setting "data_ref" to 0, thereby not adjusting
% the potential the second time this function is called. Doing this may be
% useful if a reference redox couple cannot feasibly be added during the
% desired experiment. For example, a compound may actively be decomposing
% throughout the course of the experiment - even if the reference redox
% couple were added to the electrolyte, "data_ref" may be different from
% "data". As such, it may be worthwhile to add the reference redox couple
% at the end of the experiment, and use the value of "Eref" from that final
% experiment to adjust the potential values of all the preceding
% experiments.

%% Calculate dimensionless current
phi=simcurrent(params,dataV,lib,type); % Simulate the concentration-normalized currents

%% Finding the concentration
concguess=params{1}*1000; % Convert from M to mol/m^3
% Find the initial concentration vector, as well as its upper and lower bounds
[conc,conclb,concub]=voltammogram_conc_fitting_protocol(dataInet,dataIstd,phi,ninitguess,concguess,lib); % Outputs in mol/m^3, or mmol/L
if type==1 % Assign the appropriate potential metric to Estep, whether a conventional CV
           % or CSW voltammogram is being studied.
    Estep=dataV;
else
    Estep=dataEstep;
end
inet0=phi*conc'; % Calculate inet0
BICall=biccompx(conc,phi,dataInet,dataIstd,lib,data,type); % Calculate the BIC for all compounds in solution
                                                           % (this is Model 1 from Figure 2 in the main text). 
lib0=lib; % Store library
conc0=conc; % Store concentration vectors
phi0=phi; % Store dimensionless difference currents
% The loop below loops through all the compounds in the library and
% elimates the compound being considered from contention. The resulting fit
% is then compared to that of when all compounds are considered, and the
% probability of a compound being present in the electrolyte of interest is
% calculated from this comparison.
for i=1:size(lib0,2) % Loop through the library
    lib=lib0; % Reset library, concentration vector, and conc-normalized currents from previous iteration
    conc=conc0;
    phi=phi0;
    conc(i)=[]; % Eliminate entries from the compound being examined
    lib(i)=[];
    phi(:,i)=[];
    % Fit a new concentration vector (no need for bounds)
    conc=voltammogram_conc_fitting_protocol(dataInet,dataIstd,phi,ninitguess,conc,lib);
    BIC(i)=biccompx(conc,phi,dataInet,dataIstd,lib,data,type); % Calculate the BIC for the compound set excluding the compound of interest
    BICmin=min([BIC(i),BICall]); % Find the min BIC
    BFx=exp((BICmin-BICall)/2); % Bayes Factor for compound "X" being present in the electrolyte
    BFnotx=exp((BICmin-BIC(i))/2); % Bayes Factor for compound "X" being absent in the electrolyte
    P(i)=BFx.*prior(i)./(BFx.*prior(i)+BFnotx.*(1-prior(i))); % Calculate the probability
end

%% Calculate probability that the peak isn't noise
% Estimate mean and standard deviation for concentrations for each compound
% for estimation.
%
% This section of the code makes use of concentration units to generate a
% concentration detection limit (based on background noise), taking the form of a normal distribution
% about the mean concentration. This is then compared with the estimated
% concentration of the compound being examined. If the estimated
% concentration is far from the mean + n*std deviations (n an O(1)
% integer), it is assumed that the signal attributed to the compound being
% studied has a high probability of not arising from background
% processes. However, if the estimated concentration is on the same order
% as that the mean + n*std deviations (n holding the same meaning), then it
% is probable that the signal arises from background processes - this must
% be considered when final probabilities are calculated.
for i=1:size(lib0,2)
    if conc0(i)==0 % Make the concentrations non-zero. If the concentration is zero, NaN will be generated for "stdcbg", obfsucating results.
        c(i)=1e-14; % An arbitrary small value larger than the computational rounding error
    else
        c(i)=conc0(i);
    end
    % Concentrations are not directly measured, but currents are. Currents
    % are linearly proportional to the concentration and the
    % concentration-normalized current - this relation is used to equate
    % currents with concentrations.
    [maxI,indmax]=max(phi0(:,i)*c(i)); % Calculate the max current for the given compound
    meancbg=0; % Calculate the mean estimated concentration, based on the 
               % scaling of the true concentration by the current. Since 
               % the background is subtracted, the mean current is expected
               % to be zero
    stdcbg=bgIstd(indmax)./maxI.*c(i); % Calculate the std deviation of the concentration
                                       % at the potential value at which the peak occursestimated conc
    %
    % Sample the probability using the CDF. The question is what is the probability that the observed
    % concentration is within +- that amount from the mean. Since there can't be
    % less than 0 concentration, bound the CDF from 0 to infinity. Since
    % the mean is zero, this means multiplying the CDF value by a factor of
    % two (or dividing by 1/2).
    %
    % "normcdf" calculates the integral from negative infinity to the point
    % of interest. First, calculate the integral from negative infinity to the
    % lower and upper bounds to yield a 1x2 vector (not neglecting to
    % divide by 1/2).
    Pbgvec=normcdf([0,meancbg+abs(c(i)-meancbg)],meancbg,stdcbg)./0.5;
    % Then, take the difference of the two vectors. The effect of negative
    % infinity is canceled out, yielding the CDF between the lower and
    % upper bounds. This CDF value is the probability that the observed
    % concentration is NOT background noise.
    %
    % The probability that noise could occur at or outside the bounds given
    % by the concentration is 1 - Pnobg, since x% of all noise observations will
    % take place within the bounds perscribed.
    Pnobg(i)=Pbgvec(2)-Pbgvec(1);
end
P=P.*Pnobg; % Assume the two probabilities are independent.

%% Calculate final different current estimated from the optimal current vector
conc=conc0; % Reset the concentrations, concentration-normalized currents, and the library
phi=phi0;
lib=lib0;
j=0; % Counter to count for elements that have already been deleted
for i=1:length(conc0)
    if P(i)<0.5
        conc(i-j)=[];
        conclb(i-j)=[];
        concub(i-j)=[];
        lib(i-j)=[];
        phi(:,i-j)=[];
        j=j+1;
    end
end
[concfinal,conclbfinal,concubfinal]=voltammogram_conc_fitting_protocol(dataInet,dataIstd,phi,ninitguess,conc,lib); % Find the final concentrations and bounds
inet=phi*concfinal';

%% Okay, so this works. How do we account for the variability in the probability? Account for the variability in parameters?
%% Stop timer
simtime=toc;
end

function BIC=biccompx(conc,phi,dataInet,dataIstd,lib,data,type)
%% Header
% This function calculates the BIC of a compound being present in an electrolyte, given
% the estimated concentration, the concentration-normalized current, the 
% experimental data, the standard deviation, library, and
% the original input data. All inputs are fed from
% "conc_compound_find.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2020-07-20
% Version: 1.0
%
%% Inputs:
%
% conc -     vector of estimated concentrations from the fitting routine, in
%            units of mol/m^3.
%
% phi -      concentration-normalized currents, simulated from "simcurrent.m".
%            In units of A * m^3 / mol.
%
% dataInet - the processed difference current, in units of A.
%
% dataIstd - the experimental standard deviation, accounting for deviations
%            in the background-subtracted current, in units of A.
%
% lib -      library used to identify compounds (regardless of whether the
%            compound of interest is truncated). Read the header of
%            "conc_compound_find.m" for further details.
%
% data -     the original dataset. Used to estimate the number of data points
%            used in the calculation of the BIC. Read "conc_compound_find.m"
%            for more information on formatting
%
% type -     determines whether the dataset is a conventional CV (type=1) or a
%            CSW voltammogram (type~=1).
%
%
%% Output:
%
% BIC - the output BIC value for the data being analyzed
%
%% Main body
inetsim=phi*conc'; % Simulated current
nparams=0; % Set parameter counter to zero
for i=1:size(lib,2) % Determine the number of parameters. The total number
                    % parameters is the sum of the number of parameters for
                    % each electron transfer mechanism in the library.
    if lib{i}(end)==1 % Reversible ET
        nparams=nparams+3;
    elseif lib{i}(end)==2 % Quasireversible ET
        nparams=nparams+5;
    end
end
Llikely=0; % Initialize the log likelihood function summation
for i=1:length(dataInet)
    Llikely=Llikely+log((1/((dataIstd(i))*sqrt(2*pi))))+(-(1/(2*dataIstd(i)^2))*(dataInet(i)-inetsim(i))^2); % Sum the log of the normal distribution. Use experimental standard deviation as an estimate.
end
if type==1 % Conventional CV
    BIC=log((size(data,2)/2)*length(inetsim)).*nparams-2*Llikely; % There are only two entries per dataset - E and I.
else % CSW voltammogram
    BIC=log((size(data,2)/4)*length(inetsim)).*nparams-2*Llikely; % There are four entries per dataset - baseline E, difference I, E, and I.
end
end

function phi=simcurrent(params,dataV,lib,type)
%% Header
% "simcurrent.m" simulates and generates concentration-normalized
% voltammograms from the input parameters, the potential waveform, the
% library being referenced, and the type of voltammogram being simulated
% (e.g., conventional CV vs. CSW voltammogram). The simulated
% concentration-normalized voltammograms are then output to be fitted to
% experimental data by adjusting the concentrations. All inputs are fed from
% "conc_compound_find.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2020-07-20
% Version: 1.0
%
%% Inputs:
%
% params - input parameters. Details on the values of these parameters are
%          found in the parent function.
%
% dataV -  the input experimental potential, used as the points of
%          simulation in "simcurrent.m".
%
% lib -    library of compounds. Details on the values of these parameters are
%          found in the parent function.
%
% type -   determines whether the dataset is a conventional CV (type=1) or a
%          CSW voltammogram (type~=1).
%
%% Output:
%
% phi - vector of simulated concentration-normalized currents (or difference currents
%       if a CSW voltammogram is being analyzed). Rows are the same potential
%       point, while columns represent a compound in the library. As such, "phi"
%       should be a m x nlib matrix, where there are "m" voltammetric points
%       simulated and "nlib" compounds in the library.
%
%% Main body
%% Extract parameters
concguess=params{1}; % mol/L
re=params{2}; % mm
dE=params{3}; % mV
dE=dE./1000; % mV to V
Esw=params{4}; % mV
Esw=Esw./1000; % mV to V
tpulse=params{5}; % ms
tpulse=tpulse./1000; % ms to s
T=params{6}; % K
ptspersec=params{7}; % Unitless
nuCSWVeff=2*dE/tpulse; % Effective scan rate for a CSW voltammogram
nu=params{8}; % mV/s
if type==1 % Conventional CV
    dnu=nu/1000/ptspersec; % delta-Volt per time point
else % CSW voltammogram
    dnu=nuCSWVeff/ptspersec;
end

parfor q=1:size(lib,2)
    % Extract information from the library
    lib_loop=lib{1,q};
    if lib_loop(end)==1 % Reversible ET mechanism
        Ef=lib_loop(1); % Formal potential (V vs. Ref)
        Do=lib_loop(2); % Diffusion coefficient of oxidized species (cm^2/s)
        Dr=lib_loop(3); % Diffusion coefficient of reduced species (cm^2/s)
        ko=1; % Filler value of the heterogeneous rate constant such that the code doesn't crash
        alpha=1; % Filler value of the transfer coefficient such that the code doesn't crash
    elseif lib_loop(end)==2 % Quasireversible ET mechanism
        Ef=lib_loop(1);
        Do=lib_loop(2);
        Dr=lib_loop(3);
        ko=lib_loop(4); % Heterogeneous rate constant (cm/s)
        alpha=lib_loop(5); % Transfer coefficient (unitless)
    end
    % Initialize parameters
    params=zeros(1,15);
    params(1)=concguess(q);
    params(2)=Dr;
    params(3)=Do;
    params(4)=re;
    params(5)=dE;
    params(6)=Esw;
    params(7)=tpulse;
    params(8)=T;
    params(9)=Ef;
    params(10)=alpha;
    params(11)=ko;
    params(12)=dnu;
    params(13)=nu;
    if dataV(floor(end/2))>dataV(1) % If the initial sweep is oxidative (i.e., reduced species initially present in the bulk)
        [~,~,~,~,phi(:,q)]=Voltammetry_Simulator_MacroDisk_2019_11_07(lib_loop(end),1,type,dataV,params); % Output the concentration-normalized current
    elseif dataV(floor(end/2))<dataV(1) % If the initial sweep is reductive (i.e., oxidized species initially present in the bulk)
        [~,~,~,~,phi(:,q)]=Voltammetry_Simulator_MacroDisk_2019_11_07(lib_loop(end),2,type,dataV,params);
    end
end
end

function [dataV,dataInet,Idatastd,dataEstep,bgIstd,Eref]=data_processing_conc_find_CSWV(data,bsubtract,data_ref,data_bgref,paramsin,nlib)
%% Header
% This function pre-processes raw experimental data by first 1) iR
% correcting it, 2) referencing the potentials to that of ferrocene (or any desired redox event), 
% and 3) normalizing the data so the the concentration of all the trials is the
% same. Once this adjustment is performed, then the data is averaged, the
% background current subtracted, and the standard deviation recorded for
% later in the code. Note that this code is for CSW voltammograms only;
% a separate code will pre-process the CV data. All inputs come from
% "conc_compound_find.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2020-07-20
% Version: 1.0
%
%% Inputs:
%
% data -       raw experimental data. See the parent function for details
%              on the formatting of this data.
%
% bsubtract -  the raw background current. See the parent function for details
%              on the formatting of this data.
%
% data_ref -   raw experimental data, including the reference compound, for
%              potential referencing (e.g., ferrocene) of the data. See the parent function for details
%              on the formatting of this data.
%
% data_bgref - raw experimental data, including the reference compound, for
%              potential referencing (e.g., ferrocene) of the background current. 
%              See the parent function for details on the formatting of this data. 
%
% paramsin -   the input parameters from the parent function, which contains 
%              more details on the formatting of these parameters.
%
% nlib -       the number of entries in the library. See the parent function for
%              details on the formatting of these parameters.
%
%
%% Outputs:
%
% dataV -     the raw potential values (i.e., the values of the forward and
%             reverse pulses), adjusted via iR compensation and referencing vs. the
%             redox event of interest. In V vs. ref.
%
% dataInet -  the processed difference current, normalized by the
%             concentration. In units of A.
%
% Idatastd -  the standard deviation in the difference current at a given
%             potential point over all trials. The effects of the variance in the
%             background current are also accounted for by evaluating the covariance
%             matrix of the two currents. In units of A.
%
% dataEstep - the baseline potential, adusted via iR compensation and
%             referencing vs. the ferrocene redox event. In V vs. ref.
%
% bgIstd -    the mean standard deviation of the background current. Used to
%             determine whether an observed peak occurs from the presence of an
%             electroactive compound or from background processes.
%
% Eref -      the redox potential of the reference redox couple with respect to
%             the experimental potential (V vs. Ref). Another way to look at it
%             is the amount the potentials are adjusted by.


%% Initialization
% Set a bound so that the index does not exceed the number of array
% elements when accounting for the variance in the background data.
% This bound needs to be set in place because the number of background
% trials may not be the same as the number of trials in which the only
% active species in the electrolyte is the reference redox couple (e.g.,
% ferrocene).
%
% Voltammograms containing only reference compounds are used to estimate
% the standard deviation of a single dataset.  This way, electrolytes that
% transfrom over the course of an experiment may be quantified in a single trial.
minindex=min([size(bsubtract,2)/4,size(data_bgref,2)/4]);

% Extract the dataInet, dataV, dataEstep, and dataI vectors
for i=1:size(data,2)/4
    dataInet(:,i)=data(:,4.*i-2)./1e6; % Convert from mu_A to A
    dataV(:,i)=data(:,4.*i-1);
    dataEstep(:,i)=data(:,4.*i-3);
    dataI(:,i)=data(:,4*i)./1000; % Convert from mA to A
end

% Get rid of NaN - only the first half of entries for Estep
% and Inet are significant
dataInet=dataInet(1:end/2,:);
dataEstep=dataEstep(1:end/2,:);

% Extract resistances and compensations
R=paramsin{9}; % Resistance of electrolyte being probed to identify compounds
comp=paramsin{10}; % iR compensation of electrolyte being probed to identify compounds
Rbg=paramsin{11}; % Resistance of electrolyte used to determine background signal (only supporting salt)
compbg=paramsin{12}; % iR compensation of electrolyte used to determine background signal (only supporting salt)
Rref=paramsin{13}; % Resistance of electrolyte being probed to identify compounds, including the reference redox couple
compref=paramsin{14}; % iR compensation of electrolyte being probed to identify compounds, including the reference redox couple
Rbgref=paramsin{15}; % Resistance of electrolyte with only the reference redox couple
compbgref=paramsin{16}; % iR compensation of electrolyte with only the reference redox couple
concref=paramsin{17}; % Vector of concentrations for the electrolyte with only the reference redox couple (used to normalize and estimate the standard deviation)

%% IR Correction on Raw Data
[dataEstep,dataV]=iR_correction(dataV,dataI,R,comp); % iR correct the data. There is only one dataset, so no need to loop

%% Potential Adjustment on Raw Data
% In non-aqueous solutions, ferrocene or cobaltocene can be
% used as a reference. If cobaltocene is used, this, in turn, can be
% referenced to ferrocene so that all potentials are referenced to ferrocene.
[dataEstep,dataV,Eref]=reference_correction(dataEstep,dataInet,dataV,data_ref,Rref,compref,1,nlib,data_ref);

%% Background Current Processing
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
for i=1:size(bgP,2) % There are multiple trials for the background current, so their values must be averaged
    [bgP(:,i),bgPraw(:,i)]=iR_correction(bgPraw(:,i),bgIraw(:,i),Rbg(i),compbg(i));
end

%% Adjust the potential to that vs. the desired reference potential
bgP=reference_correction(bgP,bgI,bgPraw,data_bgref,Rbgref,compbgref,0,nlib,data_bgref);

%% Average the background currents
bgPmean=mean(bgP,2);
bgImean=mean(bgI,2);
bgIstd=std(bgI,0,2);

%% Aligning the Background Current with the Raw Data
% The background data should be wider than any experimental dataset. These
% lines will select which indices of the background data correspond to the
% experimental data
distmin=abs(min(dataEstep)-bgPmean(1:(length(bgPmean)+1)/2)); % This is a 
% vector of distances between the minimum baseline potential of the data
% and the potential value of the background current in the forward sweep
[~,imin]=min(distmin); % Find the index that has the minimum distance
% between the baseline potential of the background data and the initial
% baseline potential of the data
imax=imin+(length(dataEstep)-1)/2; % The -1 is present because the first entry is excluded (it's already counted for in imin)

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
% Process the data with only the reference redox couple
for i=1:size(data_bgref,2)/4
    refPbg(:,i)=data_bgref(:,4*i-3);
    refIbg(:,i)=data_bgref(:,4*i-2)./1e6; % Convert from mu_A to A
    refPrawbg(:,i)=data_bgref(:,4*i-1);
    refIrawbg(:,i)=data_bgref(:,4*i)./1000; % Convert from mA to A
end

% Get rid of NaN or 0 - only the first half of entries for Estep
% and Inet are significant
refPbg=refPbg(1:end/2,:);
refIbg=refIbg(1:end/2,:);

%% Normalize the Currents by Concentration
for i=1:size(refIbg,2)
    refIbg(:,i)=refIbg(:,i)*1e-3./concref(i); % All currents are normalized to 1 mM
    refIrawbg(:,i)=refIrawbg(:,i).*1e-3./concref(i);
end

%% Background Correct and Estimate Standard Deviaton
A=[1,-1]; % The raw current is being added (+1), and the background current is being subtracted (-1)
for i=1:min(size(Ibgexp,1),size(refIbg,1)) % Loop over the minimum of the two lengths.
                                           % This is an okay first-pass
                                           % approximation, as the
                                           % background current is a weak
                                           % function of potential and does
                                           % not significatly change in
                                           % magnitude.
    V=cov(refIbg(i,1:minindex),Ibgexp(i,1:minindex)); % Create a covariance matrix between the raw current and background current
    Idatastd(i)=A*V*A'; % Calculate the covariance
    if Idatastd(i)<0
        error('Unphysical system - variance is negative')
    end
    Idatastd(i)=sqrt(Idatastd(i)); % Take the square root to get the std dev
end
Idatastd=Idatastd'; % Convert to a column vector
dataInet=dataInet-Iadjust;

%% Interpolate standard deviations
% This section estimates the standard deviation of the single experimental
% data trial by assuming that the standard deviation is linearly
% proportional to the magnitude of the current. This linear fit is then
% applied to the experimental data to estimate the standard deviation of
% the experimental data. Finally, the standard deviation can be emperically
% scaled to explore the effects of larger or smaller standard deviations on
% the output results.
P=polyfit(abs(mean(refIbg(1:length(Idatastd)),1))',Idatastd,1);
Idatastd=polyval(P,abs(dataInet));
Idatastd=1.*Idatastd; % The factor in front of Idatastd on the RHS is used
                      % to scale the estimated standard deviation as desired
end

function [dataEstep,dataVnew]=iR_correction(dataV,dataI,R,comp)
%% Header
% This function iR-corrects both the raw potential and the baseline
% potential, given the true current, the resistance, and the compensation.
% All inputs are directly fed from "data_processing_conc_find_CSWV.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2020-07-20
% Version: 1.0
%
%% Inputs:
% dataV - the raw (not baseline) potential values, in V vs. ref.
%
% dataI - the raw (not difference) currents, in A.
%
% R -     the solution resistance, in Ohms.
%
% comp -  the compensation (in a decimal form between 0 and 1). 1-comp is the uncompensated (e.g., at
%         85% ZIR, comp=0.85, 1-comp=0.15)
%
%
%% Outputs:
% dataEstep - the iR-corrected baseline potential, in V vs. ref
%
% dataVnew -  the iR-corrected raw potential, in V vs. ref
%
%% Main Body
dataVnew=dataV-dataI*(1-comp)*R; % Perform a ZIR correction on the raw potential
for i=1:length(dataV)/2
    dataEstep(i)=(dataV(2*i-1)-dataI(2*i-1)*(1-comp)*R+dataV(2*i)-dataI(2*i)*(1-comp)*R)/2; % Compensating the baseline potential
end
dataEstep=dataEstep';
end

function [dataEstep,dataVnew,Eref]=reference_correction(dataEstep_comp,dataInet_comp,dataV,data_ref,Rref,comp,bg,nlib,ref)
%% Header
% This function takes in potentials measured against some pseudoreference
% (e.g., Ag/Ag+) and adjusts the potential to some desired reference redox
% event (e.g., ferrocene) using data that has the redox event of interest
% (if any) and the reference redox event. This code can both change the
% reference potentials for the experimental data and the background data.
% All inputs are fed from "data_processing_conc_find_CSWV.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2020-07-20
% Version: 1.0
%
%
%% Inputs:
% dataEstep_comp - the baseline potential for the compound (if any; hence "comp")
%                  of interest without any reference compound. In units of V vs. ref.
%
% dataInet_comp -  the difference current for the compound (if any) of interest
%                  without any reference compound. In units of V vs. ref.
%
% dataV -          the raw potential including the compound (if any) of interest
%                  without any reference compound. In units of V vs. ref.
%
% data_ref -       the full matrix (i.e., currents and potentials) of the data containing the reference compound.
%                  This full matrix will be further processed in the code below.
%                  Directly fed from "data_processing_conc_find_CSWV.m", which in turn is
%                  fed from "conc_compound_find.m", where further details may be found.
%
% Rref -           the resistance of the solution. In units of Ohms.
%
% comp -           the compensation (in a decimal form between 0 and 1). 1-comp is the uncompensated (e.g., at
%                  85% ZIR, comp=0.85, 1-comp=0.15)
%
% bg -             an index determining whether the data whose potential is being
%                  corrected is the background current. Input 0 if the data is the background
%                  data, and input another value if the data is the experimental data (i.e.,
%                  containing the sets of compounds of interest).
%
% nlib -           number of entries in the library. Used to determine the maximum
%                  number of peaks that may be attributed to a compound.
%
% ref -            An indicator of whether to adjust the potential or not. If ref is the
%                  voltammogram of interest with the reference redox couple added, then
%                  the potential is adjusted. If ref is the scalar 0, then
%                  it is not. The format for the data, if the potential is
%                  being adjusted, is found in the parent function "conc_compound_find.m".
%
%
%% Outputs:
%
% dataEstep - the new baseline potential, referenced to the desired redox
%             event of interest. In V vs. ref
%
% dataVnew -  the new raw potential, referenced to the desired redox event
%             of interest. In V vs. ref
%
% Eref -      the redox potential of the reference redox couple with respect to
%             the experimental potential (V vs. Ref). Another way to look at it
%             is the amount the potentials are adjusted by.
%
%% Main body
% First, find the maximum potentials of the peak from the compound being
% tested
% The maximum number of peaks possible is 2 * number of compounds in the
% library. Start there.
for i=1:size(dataEstep_comp,2)
    [~,ind_comp,~,p_comp]=findpeaks(abs(dataInet_comp(:,i)),'SortStr','descend','NPeaks',2*nlib);
    elim_vec=0; % Pre-allocate a vector to eliminate peaks from contention that don't exhibit significant prominence
    for j=1:length(p_comp)
        if p_comp(j)<0.03*max(p_comp) % The threshold prominence is set to 3 % of the maximum (emperical standard)
            elim_vec(j)=1; % Mark for deletion
        else
            elim_vec(j)=0; % Mark to keep
        end
    end
    ind_comp(elim_vec==1)=[]; % Delete undesired entries
    p_comp(elim_vec==1)=[];
    npeaks=length(p_comp); % Find the number of peaks
    ncomps=npeaks/2; % There are two peaks per compound
    E0_comp(i)=mean([dataEstep_comp(ind_comp,i)]); % Find an agglomerate E0.
    % Although its physical significance is not readily obvious, it can be
    % used to calculate the value of the redox potential of the reference
    % couple, as will be demonstrated later
end

%% Extract the potentials and currents from the reference data if it exists
if max(max(ref))==0 % If no potential adjustment data is put in, make no adjustment
    warning('No adjustment will be made to another reference potential.')
    Eref=0;
    dataEstep=dataEstep_comp;
    dataVnew=dataV;
else
    for i=1:size(data_ref,2)/4 % Extract reference data
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
        [~,ind_ref,~,p_ref]=findpeaks(abs(refI(:,j)),'SortStr','descend','NPeaks',npeaks+2); % Evaluate the voltammograms for more peaks
        if bg==0 % Background current potential adjustment
            [~,ind_ref,~,~]=findpeaks(abs(refI(:,j)),'SortStr','descend','NPeaks',2); % Limit to two peaks if the background current is being corrected.
            Eref=mean([refP(ind_ref(1),j),refP(ind_ref(2),j)]);
            dataEstep(:,j)=dataEstep_comp(:,j)-Eref; % Adjust potentials
            dataVnew(:,j)=dataV(:,j)-Eref;
        else % Non-background current potential adjustment
            % Eliminate peaks with a prominence less than 1/0.03 times that of the peak
            % with the largest prominence. This is an emperical threshold
            % that can be adjusted if desired.
            elim_vec=0;
            for i=1:length(p_ref)
                if p_ref(i)<0.03*max(p_ref)
                    elim_vec(i)=1;
                else
                    elim_vec(i)=0;
                end
            end
            ind_ref(elim_vec==1)=[];
            p_ref(elim_vec==1)=[];
            
            if length(p_ref)==npeaks % This would only occur if the reference
                % redox event of interest has a convoluted peak with that
                % from the compound whose properties are being extracted
                warning('No adjustment will be made to another reference potential.')
                Eref=0;
                dataEstep=dataEstep_comp;
                dataVnew=dataV;
            elseif length(p_ref)==npeaks+2
                Eavg=mean(refP(ind_ref,j));
                Eref=((npeaks+2)*Eavg-npeaks*E0_comp)/2; % Knowledge of E0 can be compared with Eavg to find Eref
                                                            % This potential is the potential of the desired redox event
                                                            % vs. the pseudoreference.   
                dataEstep(:,j)=dataEstep_comp-Eref; % Adjust potentials
                dataVnew(:,j)=dataV-Eref;
            else
                error('Too many peaks identified')
            end
        end
    end
end
end

function [dataV,dataI,Idatastd,bgIstd,Eref]=data_processing_conc_find_CV(data,bsubtract,data_ref,data_bgref,paramsin,nlib)
%% Header
% This function pre-processes raw experimental data by first 1) iR
% correcting it, 2) referencing the potentials to that of ferrocene, and 3)
% normalizing the data so the the concentration of all the trials is the
% same. Once this adjustment is performed, then the data is averaged, the
% background current subtracted, and the standard deviation recorded for
% later in the code. Note that this code only pre-processes CV data - a
% separate code processes CSWV data. All inputs come from
% "conc_compound_find.m".
%
%% Inputs:
%
% data -       raw experimental data. See the parent function for details
%              on the formatting of this data.
%
% bsubtract -  the raw background current. See the parent function for details
%              on the formatting of this data.
%
% data_ref -   raw experimental data, including the reference compound, for
%              potential referencing (e.g., ferrocene) of the data. See the parent function for details
%              on the formatting of this data.
%
% data_bgref - raw experimental data, including the reference compound, for
%              potential referencing (e.g., ferrocene) of the background current. 
%              See the parent function for details on the formatting of this data. 
%
% paramsin -   the input parameters from the parent function, which contains 
%              more details on the formatting of these parameters.
%
% nlib -       the number of entries in the library. See the parent function for
%              details on the formatting of these parameters.
%
%
%% Outputs:
%
% dataV -     the raw potential values (i.e., the values of the forward and
%             reverse pulses), adjusted via iR compensation and referencing vs. the
%             redox event of interest. In units of V vs. ref.
%
% dataInet -  the processed difference current, normalized by the
%             concentration. In units of A.
%
% Idatastd -  the standard deviation in the difference current at a given
%             potential point over all trials. The effects of the variance in the
%             background current are also accounted for by evaluating the covariance
%             matrix of the two currents. In units of A.
%
% dataEstep - the baseline potential, adusted via iR compensation and
%             referencing vs. the ferrocene redox event. In units of V vs. ref.
%
% Eref -      the redox potential of the reference redox couple with respect to
%             the experimental potential (V vs. Ref). Another way to look at it
%             is the amount the potentials are adjusted by.


%% Extract relevant background information given the scan rate
% Select the background data that will be evaluated
    if paramsin{8}==25 % 25 mV/s
        bsubtract=bsubtract(:,1:6);
        data_bgref=data_bgref(:,1:6);
        paramsin{11}=paramsin{11}(1:3); % Background resistance
        paramsin{12}=paramsin{12}(1:3); % Background compensation
        paramsin{15}=paramsin{15}(1:3); % Background reference resistance 
        paramsin{16}=paramsin{16}(1:3); % Background reference compensation
        paramsin{17}=paramsin{17}(1:3); % Concentrations of ferrocene
        j=0; % Introduce a counter
        for i=1:size(bsubtract,1)
            if sum(isnan(bsubtract(i-j,:))==1) > 0 % If any elements are NaN
                bsubtract(i-j,:)=[]; % Eliminate them
                j=j+1; % Update the counter to consider the eliminated row
            end
        end
        j=0; % Repeat for the reference background current
        for i=1:size(data_bgref,1)
            if sum(isnan(data_bgref(i-j,:))==1) > 0
                data_bgref(i-j,:)=[];
                j=j+1;
            end
        end
    elseif paramsin{8}==50
        bsubtract=bsubtract(:,7:12);
        data_bgref=data_bgref(:,7:12);
        paramsin{11}=paramsin{11}(4:6); % Background resistance
        paramsin{12}=paramsin{12}(4:6); % Background compensation
        paramsin{15}=paramsin{15}(4:6); % Background reference resistance 
        paramsin{16}=paramsin{16}(4:6); % Background reference compensation
        paramsin{17}=paramsin{17}(4:6); % Concentrations of ferrocene
        j=0;
        for i=1:size(bsubtract,1)
            if sum(isnan(bsubtract(i-j,:))==1) > 0 % If any elements are NaN
                bsubtract(i-j,:)=[]; % Eliminate them
                j=j+1; % Update the counter to consider the eliminated row
            end
        end
        j=0; % Repeat for the reference background current
        for i=1:size(data_bgref,1)
            if sum(isnan(data_bgref(i-j,:))==1) > 0
                data_bgref(i-j,:)=[];
                j=j+1;
            end
        end
    elseif paramsin{8}==100
        bsubtract=bsubtract(:,13:18);
        data_bgref=data_bgref(:,13:18);
        paramsin{11}=paramsin{11}(7:9); % Background resistance
        paramsin{12}=paramsin{12}(7:9); % Background compensation
        paramsin{15}=paramsin{15}(7:9); % Background reference resistance 
        paramsin{16}=paramsin{16}(7:9); % Background reference compensation
        paramsin{17}=paramsin{17}(7:9); % Concentrations of ferrocene
        j=0;
        for i=1:size(bsubtract,1)
            if sum(isnan(bsubtract(i-j,:))==1) > 0 % If any elements are NaN
                bsubtract(i-j,:)=[]; % Eliminate them
                j=j+1; % Update the counter to consider the eliminated row
            end
        end
        j=0; % Repeat for the reference background current
        for i=1:size(data_bgref,1)
            if sum(isnan(data_bgref(i-j,:))==1) > 0
                data_bgref(i-j,:)=[];
                j=j+1;
            end
        end
    elseif paramsin{8}==200
        bsubtract=bsubtract(:,19:24);
        data_bgref=data_bgref(:,19:24);
        paramsin{11}=paramsin{11}(10:12); % Background resistance
        paramsin{12}=paramsin{12}(10:12); % Background compensation
        paramsin{15}=paramsin{15}(10:12); % Background reference resistance 
        paramsin{16}=paramsin{16}(10:12); % Background reference compensation
        paramsin{17}=paramsin{17}(10:12); % Concentrations of ferrocene
        j=0;
        for i=1:size(bsubtract,1)
            if sum(isnan(bsubtract(i-j,:))==1) > 0 % If any elements are NaN
                bsubtract(i-j,:)=[]; % Eliminate them
                j=j+1; % Update the counter to consider the eliminated row
            end
        end
        j=0; % Repeat for the reference background current
        for i=1:size(data_bgref,1)
            if sum(isnan(data_bgref(i-j,:))==1) > 0
                data_bgref(i-j,:)=[];
                j=j+1;
            end
        end
    elseif paramsin{8}==500
        bsubtract=bsubtract(:,25:30);
        data_bgref=data_bgref(:,25:30);
        paramsin{11}=paramsin{11}(13:15); % Background resistance
        paramsin{12}=paramsin{12}(13:15); % Background compensation
        paramsin{15}=paramsin{15}(13:15); % Background reference resistance 
        paramsin{16}=paramsin{16}(13:15); % Background reference compensation
        paramsin{17}=paramsin{17}(13:15); % Concentrations of ferrocene
        j=0;
        for i=1:size(bsubtract,1)
            if sum(isnan(bsubtract(i-j,:))==1) > 0 % If any elements are NaN
                bsubtract(i-j,:)=[]; % Eliminate them
                j=j+1; % Update the counter to consider the eliminated row
            end
        end
        j=0; % Repeat for the reference background current
        for i=1:size(data_bgref,1)
            if sum(isnan(data_bgref(i-j,:))==1) > 0
                data_bgref(i-j,:)=[];
                j=j+1;
            end
        end
    elseif paramsin{8}==1000 % For 1000 mV/s, there are only two background
                           % datasets considered (one is excluded because
                           % it's a poor dataset)
        bsubtract=bsubtract(:,31:34);
        data_bgref=data_bgref(:,31:34);
        paramsin{11}=paramsin{11}(16:17); % Background resistance
        paramsin{12}=paramsin{12}(16:17); % Background compensation
        paramsin{15}=paramsin{15}(16:17); % Background reference resistance 
        paramsin{16}=paramsin{16}(16:17); % Background reference compensation
        paramsin{17}=paramsin{17}(16:17); % Concentrations of ferrocene
        j=0;
        for i=1:size(bsubtract,1)
            if sum(isnan(bsubtract(i-j,:))==1) > 0 % If any elements are NaN
                bsubtract(i-j,:)=[]; % Eliminate them
                j=j+1; % Update the counter to consider the eliminated row
            end
        end
        j=0; % Repeat for the reference background current
        for i=1:size(data_bgref,1)
            if sum(isnan(data_bgref(i-j,:))==1) > 0
                data_bgref(i-j,:)=[];
                j=j+1;
            end
        end
    end

%% Initialization
% Set a bound so that the index does not exceed the number of array
% elements when accounting for the variance in the background data.
% This bound needs to be set in place because the number of background
% trials may not be the same as the number of trials in which the only
% active species in the electrolyte is the reference redox couple (e.g.,
% ferrocene).
%
% Voltammograms containing only reference compounds are used to estimate
% the standard deviation of a single dataset.  This way, electrolytes that
% transfrom over the course of an experiment may be quantified in a single trial.
minindex=min([size(bsubtract,2)/2,size(data_bgref,2)/2]);

% Extract the dataV and dataI vectors
for i=1:size(data,2)/2
    dataV(:,i)=data(:,2.*i-1);
    dataI(:,i)=data(:,2*i)./1000; % Convert from mA to A
end

% Extract resistances and compensations
R=paramsin{9}; % Resistance of electrolyte being probed to identify compounds
comp=paramsin{10}; % iR compensation of electrolyte being probed to identify compounds
Rbg=paramsin{11}; % Resistance of electrolyte used to determine background signal (only supporting salt)
compbg=paramsin{12}; % iR compensation of electrolyte used to determine background signal (only supporting salt)
Rref=paramsin{13}; % Resistance of electrolyte being probed to identify compounds, including the reference redox couple
compref=paramsin{14}; % iR compensation of electrolyte being probed to identify compounds, including the reference redox couple
Rbgref=paramsin{15}; % Resistance of electrolyte with only the reference redox couple
compbgref=paramsin{16}; % iR compensation of electrolyte with only the reference redox couple
concref=paramsin{17}; % Vector of concentrations for the electrolyte with only the reference redox couple (used to normalize and estimate the standard deviation)

%% IR Correction on Raw Data
for i=1:size(dataV,2)
    dataV(:,i)=dataV(:,i)-dataI(:,i)*(1-comp(i))*R(i); % No additional function needed - the correction is simple enough
end

%% Potential Adjustment on Raw Data
% In non-aqueous solutions, ferrocene or cobaltocene can be
% used as a reference. If cobaltocene is used, this, in turn, can be
% referenced to ferrocene so that all potentials are referenced to ferrocene.
[dataV,Eref]=reference_correction_CV(dataV,dataI,data_ref,Rref,compref,1,nlib,data_ref);

%% Background Current Processing
% Process background subtraction vectors
for i=1:size(bsubtract,2)/2
    bgP(:,i)=bsubtract(:,2*i-1);
    bgI(:,i)=bsubtract(:,2*i)./1000; % Convert from mA to A
end

%% iR-Correct Background Data
for i=1:size(bgP,2)
    bgP(:,i)=bgP(:,i)-bgI(:,i)*(1-compbg(i))*Rbg(i);
end

%% Adjust the potential to that vs. the desired reference potential
bgP=reference_correction_CV(bgP,bgI,data_bgref,Rbgref,compbgref,0,nlib,data_bgref);

%% Average the background currents
bgPmean=mean(bgP,2);
bgImean=mean(bgI,2);
bgIstd=std(bgI,0,2);

%% Aligning the Background Current with the Raw Data
% The background data should be wider than any experimental dataset. These
% lines will select which indices of the background data correspond to the
% experimental data. Explanation behind the thought process for this code
% is fully expounded upon in the CSWV analogue of this function.

% In the case of CV, it is better to define the forward and reverse bounds
% of the background data on their own
indbgmax=floor(length(bgPmean)/2); % Define the half-way point in the background data vector
distmin=abs(min(dataV(1:floor(length(dataV)/2)))-bgPmean(1:indbgmax)); % Find where the background
                                                                       % data and the data
                                                                       % to be examined intersect
                                                                       % in the forward sweep 
distmin1=abs(max(dataV(floor(length(dataV)/2)+1:end))-bgPmean(indbgmax+1:end)); % The same for the reverse sweep 
[~,imin]=min(distmin); % Find the index of intersection for the forward sweep
[~,imin1]=min(distmin1); % The same as the line above for the reverse sweep
imax=imin+floor(length(dataV)/2); % Find the final index from the background data for the forward sweep
imax1=imin1+floor(length(dataV)/2); % Perform the same for the reverse sweep

% Extract the background current to be adjusted
Iadjust=[bgImean(imin:imax);bgImean(indbgmax+1+imin1:indbgmax+1+imax1)];
Ibgexp=[bgI(imin:imax,:);bgI(indbgmax+1+imin1:indbgmax+1+imax1,:)];
% imax-2, rather than imax-1, is used because the
% Biologic potentiostat reaches the max baseline potential on the forward
% sweep. It does not repeat the max potential on the reverse sweep. For
% example, the baseline potential waveform is of the form 0.5 V, 0.6 V,
% 0.7 V, 0.6 V, etc...; instead of 0.5 V, 0.6 V, 0.7 V, 0.7 V, 0.6 V,...
% Sometimes, the indices don't align, resulting in an error. The "if"
% statement below addressed the cases found in this work; additional
% statements may need to be added if different datasets are being examined.
if length(Iadjust)==length(dataV)-1
    Iadjust=[bgImean(imin:imax);bgImean(indbgmax+imin1:indbgmax+1+imax1)];
    Ibgexp=[bgI(imin:imax,:);bgI(indbgmax+imin1:indbgmax+1+imax1,:)];
elseif length(Iadjust)==length(dataV)+1
    Iadjust=[bgImean(imin:imax);bgImean(indbgmax+2+imin1:indbgmax+1+imax1)];
    Ibgexp=[bgI(imin:imax,:);bgI(indbgmax+2+imin1:indbgmax+1+imax1,:)];
elseif length(Iadjust)==length(dataV)+2
    Iadjust=[bgImean(imin:imax);bgImean(indbgmax+3+imin1:indbgmax+1+imax1)];
    Ibgexp=[bgI(imin:imax,:);bgI(indbgmax+3+imin1:indbgmax+1+imax1,:)];
end

%% Standard Deviation Calculation
% Process the reference data
for i=1:size(data_bgref,2)/2
    refPbg(:,i)=data_bgref(:,2*i-1);
    refIbg(:,i)=data_bgref(:,2*i)./1000; % Convert from mA to A
end

%% Normalize the Currents by Concentration
for i=1:size(refIbg,2)
    refIbg(:,i)=refIbg(:,i)*1e-3./concref(i); % All currents of voltammograms with only the reference redox couple are normalized to 1 mM
end

A=[1,-1]; % The raw current is being added (+1), and the background current is being subtracted (-1)
for i=1:min([size(refIbg,1),size(Ibgexp,1)])
    V=cov(refIbg(38+i,1:minindex),Ibgexp(i,1:minindex)); % Create a covariance matrix between the raw current and background current
    Idatastd(i)=A*V*A'; % Calculate the covariance
    if Idatastd(i)<0
        error('Unphysical system - variance is negative')
    end
    Idatastd(i)=sqrt(Idatastd(i)); % Take the square root to get the std dev
end
Idatastd=Idatastd'; % Convert to a column vector
dataI=dataI-Iadjust;

%% Interpolate standard deviations
% Explanations and assumptions are listed in the CSWV analogue of this
% code, as the thought process is similar
P=polyfit(abs(mean(refIbg(39:length(Idatastd)+38),1))',Idatastd,1);
Idatastd=polyval(P,abs(dataI)); % Note that the standard deviation is approximately zero when the current is approximately zero.
                                % This is experimentally observed for CSWV, but is not necessarily the case for CV
Idatastd=1.*Idatastd-2.*min([0,min(Idatastd)]); % Make sure the linear fit is always positive via vertical scaling
end

function [dataVnew,Eref]=reference_correction_CV(dataV,dataI,data_ref,Rref,comp,bg,nlib,ref)
%% Header
% This function takes in potentials measured against some pseudoreference
% (e.g., Ag/Ag+) and adjusts the potential to some desired reference redox
% event (e.g., ferrocene) using data that has the redox event of interest
% (if any) and the reference redox event. This code can both change the
% reference potentials for the experimental data and the background data.
% All inputs are fed from "data_processing_conc_find_CV.m".
%
% -------------
% Author:  Alexis M. Fenton Jr.
% Date:    2020-07-20
% Version: 1.0
%
%
%% Inputs:
% dataV -    the potential for the experimental data without any reference
%            compound present. In units of V vs. ref.
%
% dataI -    the current for the experimental data. In units of A.
%
% data_ref - the full matrix (i.e., currents and potentials) of all data containing the reference compound.
%            This full matrix will be further processed in the code below.
%            Directly fed from "data_processing_library_find.m", which in turn is
%            fed from "library_entry_params_compound.m", where further details may be found.
%
% Rref -     the resistance of the solution. In units of Ohms.
%
% comp -     the compensation (in a decimal form between 0 and 1). 1-comp is the uncompensated (e.g., at
%            85% ZIR, comp=0.85, 1-comp=0.15)
%
% bg -       an index determining whether the data whose potential is being
%            corrected is the background current. Input 0 if the data is the background
%            data, and input another value if the data is the experimental data (i.e.,
%            containing the compound of interest).
%
% ref -      An indicator of whether to adjust the potential or not. If ref is the
%            voltammogram of interest with the reference redox couple added, then
%            the potential is adjusted. If ref is the scalar 0, then
%            it is not. The format for the data, if the potential is
%            being adjusted, is found in the parent function "conc_compound_find.m".
%
%% Outputs:
%
% dataVnew - the new raw potential, referenced to the desired redox event
%            of interest. In V vs. ref
%
% Eref -     the redox potential of the reference redox couple with respect to
%            the experimental potential (V vs. Ref). Another way to look at it
%            is the amount the potentials are adjusted by.
%
%% Main Body
% The process for this code is analogous to its CSWV counterpart.
% Divergences in the overall process will be highlighted here - otherwise,
% refer to the CSWV counterpart for high-level concepts.
for i=1:size(dataV,2)
    % A significant difference between CV and CSWV data is that the
    % potential mesh for CV data is much more fine (i.e., the change in
    % potential between two recorded points is much smaller) and noiser.
    % The latter poses a large problem for adjusting the potential, as each
    % peak an experimentalist would identify as a single peak is actually
    % comprised of multiple small peaks (due to noise). As such, the data
    % must be smoothed so the protocol may identify peaks in a manner
    % similar to an experimentalist.
    %
    % Further, the current response as the sweep direction shifts (either
    % from oxidative to reductive or vice versa) can cause unexpected
    % behavior in the current response, resulting in false optima. To
    % overcome this issue, the maxima from the initial sweep are found, and
    % then the maxima from the reverse sweep (the current response is
    % multiplied by negative 1) are found. Naturally, this process is for
    % an initial oxidative sweep - the sign conventions would change if the
    % initial sweep were reductive.
    indmax=floor(length(dataI)/2); % Find the index in which the sweep direction changes
    % Since there may be many local maxima, a minimum peak distance (in
    % units of indices) is established to select the absolute maximum
    [~,ind_1,~,p_1]=findpeaks(smooth(dataI(1:indmax,i)),'SortStr','descend','NPeaks',nlib,'MinPeakDistance',20); % Initial oxidative sweep
    [~,ind_2,~,p_2]=findpeaks(-smooth(dataI(indmax+1:end,i)),'SortStr','descend','NPeaks',nlib,'MinPeakDistance',20); % Final reductive sweep
    ind_comp=[ind_1;ind_2+indmax]; % Concatonate the findings
    p_comp=[p_1;p_2];
    
    % Similar process as the CSWV code
    elim_vec=0;
    for j=1:length(p_comp)
        if p_comp(j)<0.01*max(p_comp) % The threshold is 1 % instead of 3 %.
                                      % Of course, this is an emperical
                                      % threshold. This threshold is set
                                      % lower because prominence of CV
                                      % peaks are generally lower than that
                                      % of CSWV peaks (as the baseline
                                      % doesn't reach zero as rapidly in
                                      % the case of CV).
            elim_vec(j)=1;
        else
            elim_vec(j)=0;
        end
    end
    ind_comp(elim_vec==1)=[];
    p_comp(elim_vec==1)=[];
    npeaks=length(p_comp);
    ncomps=npeaks/2;
    E0_comp(i)=mean([dataV(ind_comp,i)]);
end

%% Extract the potentials and currents from the reference data if it exists
if max(max(ref))==0 % If no potential adjustment data is put in, make no adjustment
    warning('No adjustment will be made to another reference potential.')
    Eref=0;
    dataVnew=dataV;
else
    % Process reference potential data
    for i=1:size(data_ref,2)/2
        refP(:,i)=data_ref(:,2*i-1);
        refI(:,i)=data_ref(:,2*i)./1000; % Convert from mA to A
    end
    
    %% iR-Correct the reference data
    for i=1:size(refP,2)
        refP(:,i)=refP(:,i)-refI(:,i)*(1-comp(i))*Rref(i);
    end
    
    %% Integrate the reference data
    % Loop over all reference datasets here
    for j=1:size(refP,2)
        % Similar concept as above (treat each sweep separately,
        % concatonate the results)
        indrefmax=floor(length(refI)/2);
        [~,ind_ref1,~,p_ref1]=findpeaks(smooth(refI(1:indrefmax,i)),'SortStr','descend','NPeaks',npeaks/2+1,'MinPeakDistance',20);
        [~,ind_ref2,~,p_ref2]=findpeaks(-smooth(refI(indrefmax+1:end,i)),'SortStr','descend','NPeaks',npeaks/2+1,'MinPeakDistance',20);
        ind_ref=[ind_ref1;ind_ref2+indrefmax];
        p_ref=[p_ref1;p_ref2];
        if bg==0 % Background current potential adjustment
            indrefmax=floor(length(refI)/2);
            [~,ind_ref1]=findpeaks(smooth(refI(1:indrefmax,i)),'SortStr','descend','NPeaks',1,'MinPeakDistance',20);
            [~,ind_ref2]=findpeaks(-smooth(refI(indrefmax+1:end,i)),'SortStr','descend','NPeaks',1,'MinPeakDistance',20);
            ind_ref=[ind_ref1;ind_ref2+indrefmax];
            Eref=mean([refP(ind_ref(1),j),refP(ind_ref(2),j)]);
            dataVnew(:,j)=dataV(:,j)-Eref; % Adjust potential
        else
            elim_vec=0;
            for i=1:length(p_ref)
                if p_ref(i)<0.01*max(p_ref)
                    elim_vec(i)=1;
                else
                    elim_vec(i)=0;
                end
            end
            ind_ref(elim_vec==1)=[];
            p_ref(elim_vec==1)=[];            
            if length(p_ref)==npeaks
                warning('No adjustment will be made to another reference potential.')
                Eref=0;
                dataVnew=dataV;
            elseif length(p_ref)==npeaks+2
                Eavg=mean(refP(ind_ref,j));
                Eref=((npeaks+2)*Eavg-npeaks*E0_comp(j))/2;
                dataVnew(:,j)=dataV(:,j)-Eref;
            else
                error('Too many peaks identified')
            end
        end
    end
end
end