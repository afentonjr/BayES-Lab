function [t,E]=E_find_interpolate(dataV,sub,tau)
%% Header
% This function takes the input potential waveform and expands it to have
% the number of elements specified by "sub" by filling in the gaps between
% the time points in the vector of potential values. This is useful for simulating
% voltammograms when the time mesh must be very fine.
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Inputs:
%
% dataV - the input potential waveform, in V vs. ref.
%
% sub -   the number of substeps (i.e., the number of time, or potential,
%         steps).
%
% tau -   the time duration of each potential, in s. In the case of CV, this
%         results in a staircase waveform. In a CSWV, this results in the
%         square-wave potential waveform.
%
%
%% Outputs:
%
% t - the expanded time vector with the number of elements specified by
%     "sub". In units of s
%
% E - the expanded potential vector with the number of elements specified
%     by "sub". In units of V vs. ref.
%
%
%% Main Body
% Make the voltage data the input potential
lengthdata=length(dataV); % Determine the number of data points in the input data
subvec=linspace(1,sub,sub); % Create a vector from 1 to the number of desired subdivisions
npperdata=sub./lengthdata; % Define the number of points at a given potential
% This loop will generate the potential/time relation assuming that all
% potentials as a function of time are equal to the reported potential at
% the end of the potential hold.
for i=1:sub
    term=ceil(subvec(i)./npperdata); % Define the index that helps denote when the potential will change
    if term>lengthdata % This bounds what indices of dataV are accessed so no error message is displayed
        % (otherwise, index would exceed array bounds).
        E(i)=dataV(term-1);
    else
        E(i)=dataV(term);
    end
end
t=tau.*lengthdata.*subvec./sub; % Make t a vector also. Assumes equal times between each potential hold
E=E';
t=t';
end