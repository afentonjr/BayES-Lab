% BayES' Lab Example Script
% Library development for PT: this script provides an example on how to
% create a library entry for a particular compound (in this case,
% unsubstituted phenothiazine).

% This file is a part of the BayES' Lab toolbox
%
% Official GitHub: """""""""""
%
% BayES' Lab: A MATLAB framework used to identify electroactive components
% from a voltammogram using a library of electrochemical-transport
% descriptors.

% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0

% Clear workspace
clear all
clc
close all

% Add to path library development directory
addpath('Library Development');

% Load files
dataPT=load('Library Development\PT Inputs and Previous Results\2020-06-30-PT-Data-for-Param-Estimation.mat');
params_inPT=load('Library Development\PT Inputs and Previous Results\2020-06-30-Input-Params-for-PT.mat');
params_guessPT=load('Library Development\PT Inputs and Previous Results\2020-06-30-Guess-Params-for-PT.mat');
conc_inPT=load('Library Development\PT Inputs and Previous Results\2020-06-30-Input-Concentrations-for-PT.mat');
bsubtract=load('Library Development\PT Inputs and Previous Results\2020-06-30-Background-Data-for-Background-Subtraction.mat');
dataPTref=load('Library Development\PT Inputs and Previous Results\2020-06-30-PT-Data-with-Ferrocene-for-Param-Estimation-Potential-Reference.mat');
dataFc_trunc=load('Library Development\PT Inputs and Previous Results\2020-06-30-Background-Data-for-Background-Subtraction-Potential-Reference.mat');
priorPT=load('Library Development\PT Inputs and Previous Results\2020-07-21-Priors-for-PT-Parameter-and-ETM-Estimation.mat');

% Expand structure
dataPT=dataPT.dataPT;
params_inPT=params_inPT.params_inPT;
params_guessPT=params_guessPT.params_guessPT;
conc_inPT=conc_inPT.conc_inPT;
bsubtract=bsubtract.bsubtract;
dataPTref=dataPTref.dataPTref;
dataFc_trunc=dataFc_trunc.dataFc_trunc;
priorPT=priorPT.priorPT;

% Execute function
params_inPT{6}=200;
[paramsfoundPT,paramslbPT,paramsubPT,ETMPT,PETMPT,dataEstepPT,dataInetPT,inetsimPT,ErefPT,simtimePT]=library_entry_params_compound(dataPT,params_inPT,params_guessPT,conc_inPT,0,bsubtract,dataPTref,dataFc_trunc,priorPT);

% Plot results of fit
% Plot the post-processed data and the best overall fit
plot(dataEstepPT,dataInetPT,'ob','LineWidth',2)
hold on
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
xlabel('Baseline Potential (V vs. FeCp_2^{0/+})')
ylabel('Difference Current (A)')
if ETMPT==1
    plot(dataEstepPT,inetsimPT(:,2),'-k','LineWidth',2)
    plot(dataEstepPT,inetsimPT(:,1),'-r','LineWidth',2)
    legend('Data','Quasi (Less Likely)','Rev (More Likely)','Location','Best')
elseif ETMPT==2
    plot(dataEstepPT,inetsimPT(:,1),'-r','LineWidth',2)
    plot(dataEstepPT,inetsimPT(:,2),'-k','LineWidth',2)
    legend('Data','Rev (Less Likely)','Quasi (More Likely)','Location','Best')
end