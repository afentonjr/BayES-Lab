function [conctemp,lbfind,ubfind]=voltammogram_conc_fitting_protocol(dataInet,dataIstd,phi,ninitguess,concguess,lib)
%% Header
% This function takes in the post-processed current, the estimated standard
% deviation, the simulated concentration-normalized current, the number of
% initial guesses, the initial concentration guesses, as well as the
% library being referenced, to estimate the concentrations of each of the
% compounds of each of the library entries, as well as its upper and lower
% bounds, using linear weighted least squares fitting with a constraint (i.e.,
% that concentrations are non-negative).
% This function is called upon by "conc_compound_find.m"; all inputs come
% from this parent function as well.
%
% -------------
% Authors: Alexis M. Fenton Jr., Fikile R. Brushett
% Date:    2021-04-13
% Version: 1.0
%
%% Inputs:
%
% dataInet -   the post-processed current, in units of A. See the parent code
%              for further output information.
%
% dataIstd -   the estimated standard deviation in the current for the
%              experimental data of interest, in units of A. See the parent code
%              for further output information.
%
% phi -        the concentration-normalized current generated by "simcurrent.m".
%              See this function for more information on phi. In units of A * m^3 / mol
%
% ninitguess - the number of initial guesses this protocol will repeat.
%
% concguess -  the initial guess for the concentrations, in mol/m^3 or mmol/L. See more
%              information in the parent code.
%
% lib -        the library being referenced for this study. See more information
%              in the parent code.
%
%
%% Outputs:
%
% conctemp -   the vector of estimated best-fit concentrations for all
%              compounds in the library. In mol/m^3 or mmol/L.
%
% lbfind -     the vector of the estimated lower bound for the best-fit concentrations for all
%              compounds in the library according to the metric of choice. In mol/m^3 or mmol/L.
%
% ubfind -     the vector of the estimated upper bound for the best-fit concentrations for all
%              compounds in the library according to the metric of choice. In mol/m^3 or mmol/L.
%
%% Main body
% Define lower and upper bounds for concentrations
lowerb=zeros(1,size(lib,2)); % Cannot have non-negative concentrations
upperb=100000.*ones(1,size(lib,2)); % Cannot feasibly have a solution of 100 M
cguessend=concguess;
% Pre-allocate vectors that will contain the estimates of the best-fit
% concentrations
cc=zeros(1,size(lib,2));
dd=zeros(1,size(lib,2));
for i=1:ninitguess+1
    % Use a rng to generate an initial concentration guess using a log
    % scale.
    randco=ceil(11*rand);
    if randco==1
        concguess(i,:)=1e-6+rand(size(lowerb)).*(1e-5-1e-6);
    elseif randco==2
        concguess(i,:)=1e-5+rand(size(lowerb)).*(1e-4-1e-5);
    elseif randco==3
        concguess(i,:)=1e-4+rand(size(lowerb)).*(1e-3-1e-4);
    elseif randco==4
        concguess(i,:)=1e-3+rand(size(lowerb)).*(1e-2-1e-3);
    elseif randco==5
        concguess(i,:)=1e-2+rand(size(lowerb)).*(1e-1-1e-2);
    elseif randco==6
        concguess(i,:)=1e-1+rand(size(lowerb)).*(1e0-1e-1);
    elseif randco==7
        concguess(i,:)=1e0+rand(size(lowerb)).*(1e1-1e0);
    elseif randco==8
        concguess(i,:)=1e1+rand(size(lowerb)).*(1e2-1e1);
    elseif randco==9
        concguess(i,:)=1e2+rand(size(lowerb)).*(1e3-1e2);
    elseif randco==10
        concguess(i,:)=1e3+rand(size(lowerb)).*(1e4-1e3);
    elseif randco==11
        concguess(i,:)=1e4+rand(size(lowerb)).*(1e5-1e4);    
    end % randomize initial X for solver
    if i==ninitguess+1
        concguess(i,:)=cguessend; % If the last initial guess is being pulled, use the input initial guess
    end
end
for i=1:ninitguess+1
    cguess=concguess(i,:);
    options=optimoptions(@fmincon,'Algorithm','interior-point','FiniteDifferenceType','central','Display','none','StepTolerance',1e-12,'TolFun',1e-12,...
        'MaxFunctionEvaluations',500,'SubproblemAlgorithm','cg'); % Options for fmincon - conjugate gradient (cg)
    [cc(i,:),fval1(i)]=fmincon(@(c)obj_func_min_conc_fit_CSWV(c,phi,dataInet,dataIstd),...
        cguess,[],[],[],[],lowerb,...
        upperb,[],options);
    options1=optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','central','Display','none','StepTolerance',1e-12,'TolFun',1e-12,...
        'MaxFunctionEvaluations',500); % Options for fmincon - sequential quadratic programming (sqp)
    [dd(i,:),fval2(i)]=fmincon(@(c)obj_func_min_conc_fit_CSWV(c,phi,dataInet,dataIstd),...
        cguess,[],[],[],[],lowerb,...
        upperb,[],options1);
    subcycleperc=i/(ninitguess+1) % Counter to track progress
end
[min1,ind1]=min(fval1); % Find the minimum of all initial guesses for cg
[min2,ind2]=min(fval2); % Find the minimum of all initial guesses for sqp
if min1>min2 % Find the overall minimum between cg and sqp
    conctemp=dd(ind2,:);
    options2=optimoptions(@fmincon,'Algorithm','sqp','FiniteDifferenceType','central','Display','none','StepTolerance',1e-12,'TolFun',1e-12,...
            'MaxFunctionEvaluations',500); % For the determination of the bounds, use the set of options that found the best fit (to save time)
elseif min1<=min2
    conctemp=cc(ind1,:);
    options2=optimoptions(@fmincon,'Algorithm','interior-point','FiniteDifferenceType','central','Display','none','StepTolerance',1e-12,'TolFun',1e-12,...
            'MaxFunctionEvaluations',500,'SubproblemAlgorithm','cg');
end
parfor q=1:length(conctemp) % Parallelize over the lower bounds first
    lb=lowerb(q); % The lower bound is the same as for the parameter estimation protocol
    ub=0.999*conctemp(q); % The upper bound is 99.9 % the estimate of the parameters from the previous fitting protocol. This extrema is typically not reached
    lbfind(q)=fmincon(@(c)obj_func_conc_bounds(c,phi,dataInet,dataIstd,conctemp,min([min1,min2]),q),...
        conctemp(q),[],[],[],[],lb,ub,[],options2);
end
parfor q=1:length(conctemp) % Then parallelize over the upper bounds
    lb=conctemp(q); % The lower bound here is the estimate of the parameters from the previous fitting protocol
    ub=upperb(q); % Keep the original bound
    ubfind(q)=fmincon(@(c)obj_func_conc_bounds(c,phi,dataInet,dataIstd,conctemp,min([min1,min2]),q),...
        conctemp(q),[],[],[],[],lb,ub,[],options2);
end
end

function f=obj_func_conc_bounds(c,phi,dataInet,dataIstd,conc,fval,q)
%% Header
% Objective function to find the best estimate of the lower and upper parameter bounds, 
% given the experimental data and how well the best estimate of the parameters
% fit the experimental data. All data imported from
% "voltammogram_conc_fitting_protocol.m".
%
%
%% Inputs:
%
% c -        the lower or upper bound to be estimated. Units of mol/m^3.
%
% phi -      the concentration-normalized current. More information may be found
%            in "simcurrent.m". Units of A * m^3 / mol.
%
% dataInet - the processed experimental current. Units of A.
%
% dataIstd - the standard deviation of the experimental current, including
%            the standard deviation from the background current. Units of A.
%
% conc -     the vector of estimated best-fit concentrations from
%            "voltammogram_conc_fitting_protocol.m". Units of mol/m^3
%
% fval -     the value of the minimized argument in the exponent. Unitless
%
% q -        the index which is used to determine which parameter the lower or
%            upper bound is being estimated for.
%
%
%% Outputs:
%
% f - the objective function to be minimized (in this case, the difference
%     of the likelihood given an estimate of the bound and half of the maximum
%     likelihood).
%
%% Main body
conc(q)=c; % Set the value of the concentration vector to be that of the
           % parameter being fitted. The concentrations of all other
           % compounds are held at their optimal value
% if fittype==0 % Used for maximizing the likelihood function, for a consistent BIC calculation
    ivec=(dataInet-phi*conc').^2./dataIstd.^2; % Calculate the argument of the normal distribution
    isum=sum(ivec);
    f=(isum-(fval-log(0.5))).^2; % The 0.5 is indicative that the likelihood function is halved.
% else % Simple least squares fitting. Seems to be better for concentration estimation
%     ivec=dataInet-phi*conc';
%     isum=norm(ivec);
%     f=(isum-0.5*fval).^2; % Half the sum of the errors.
% end
% The square is minimized, which is equivalent to zeroing the objective function
end