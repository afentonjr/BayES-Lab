function f=obj_func_min_conc_fit_CSWV(c,phi,dataInet,dataIstd)
%% Header
% Objective function to find the best estimate of the vector of best-fit concentrations
% for all the compounds in the library, given the experimental data and
% concentration-normalized difference currents. All data imported from
% "voltammogram_conc_fitting_protocol.m".
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Inputs:
%
% c -        the vector of best-fit concentrations to be estimated. Units of
%            mol/m^3.
%
% phi -      the concentration-normalized current. More information may be found
%            in "simcurrent.m". Units of A * m^3 / mol.
%
% dataInet - the processed experimental current. Units of A.
%
% dataIstd - the standard deviation of the experimental current, including
%            the standard deviation from the background current. Units of A.
%
%% Outputs:
%
% f - the objective function to be minimized (in this case, the difference
%     of the likelihood given an estimate of the bound and half of the maximum
%     likelihood).
%
%% Main body
% if fittype==0 % Used for maximizing the likelihood function, for a consistent BIC calculation
    ivec=(dataInet-phi*c').^2./dataIstd.^2; % The argument of the exponential term of the likelihood function
% %                                         (i.e., a multivariate normal distribution) is minmiized to
% %                                         maximize the likelihood function.   
    f=sum(ivec);
% else % Simple least squares fitting. Seems to be better for concentration estimation
%     ivec=dataInet-phi*c'; % Minimize the error
%     f=norm(ivec);
% end
end