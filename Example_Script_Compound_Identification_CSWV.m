% BayES' Lab Example Script
% Compound identification: this script provides an example on how to label
% compounds from an example voltammogram when provided with a library using
% cyclic square wave voltammetry (CSWV).

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
addpath('Compound Identification');

% Load files
PTMPTdatafirstelyte=load('Compound Identification\2021-03-23-PT-MPT-CSWV-Data-All-Trials-First-Electrolyte-Adjusted-Potential.mat');
params_indata=load('Compound Identification\2020-07-21-Input-Params-CSWV-Compound-Identify.mat');
library=load('Compound Identification\2021-04-13-Phenothiazine-Library.mat');
bsubtract=load('Compound Identification\2020-06-30-Background-Data-for-Background-Subtraction.mat');
dataFc_trunc=load('Compound Identification\2020-06-29-Truncated-Ferrocene-Data.mat');
prior=load('Compound Identification\2020-07-21-Prior-List.mat');

% Expand structure
PTMPTdatafirstelyte=PTMPTdatafirstelyte.PTMPTdatafirstelyte;
params_indata=params_indata.params_indata;
library=library.library;
bsubtract=bsubtract.bsubtract;
dataFc_trunc=dataFc_trunc.dataFc_trunc;
prior=prior.prior;

% Execute function
[P,conc0,concfinal,conclbfinal,concubfinal,Estep,inet0,inet,dataInet,Eref,simtime]=conc_compound_find(PTMPTdatafirstelyte(:,1:4),params_indata,99,library,bsubtract,0,dataFc_trunc,prior,0);

% Plot probabilities as a bar graph
xlab=categorical({'PT','MPT','EPT','iPrPT','PhPT'});
xlab=reordercats(xlab,{'PT','MPT','EPT','iPrPT','PhPT'});
bar(xlab,P,'LineWidth',2)
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
ylim([0 1.2])
ylabel('Posterior Probability')

% Plot the post-processed data and the best overall fit
figure()
plot(Estep,dataInet,'ob','LineWidth',2)
hold on
plot(Estep,inet,'-r','LineWidth',2)
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
xlabel('Baseline Potential (V vs. FeCp_2^{0/+})')
ylabel('Difference Current (A)')
legend('Data','Best Fit','Location','Best')

% Plot estimated concentrations as a bar graph
xlab=xlab(P>0.5); % Determine which abcissa labels should be included
figure()
bar(xlab,concfinal,'LineWidth',2)
set(gca,'FontName','Arial','FontSize',15,'LineWidth',2)
ylabel('Estimated Concentration (mM)')