function [pars covpars SE CI Zscores pvalues gval exitflag hess]=aft(time,Z,status,distr,init,minimizer,options,printresults)
%aft(surv_time(randset1), [ones(round(n*0.8),1),agt20],surv_stat(randset1),'weibull');

% The “aft” function fits models of the form:
%
% Y=log(T)=g0+g1*Z1+g2*Z2+...+sigma*epsilon
%
% where usually T is a time to event variable and g0, g1, ... and sigma are to
% be estimated. Since T is a time to event variable censoring might be involved. 
% The “aft” function deals with possibly right and/or left censored data. With "sigma" we denote the scale parameter, 
% and the regression coefficients are denoted by vector g=[g0 g1 g2...]. The design (covariate)
% matrix is Z=[Z1 Z2...], where Z1,Z2,... are the corresponding covariate column vectors, 
% (each covariate is one column). The user must include an additioanl column of ones if an intercept is to
% included. The distribution of "epsilon" defines the distribution of T. The user can specify
% the distribution with the following available options: 
% Exponential, Weibull, Log-normal, Log-logistic, Generalized Gamma. 
%
% The “aft” routine is supposed to be a MATLAB alternative to proc lifereg of SAS, or survreg of R. 
% However the “aft” has less options.
%
% For the supported distributions of “aft” and the details of the underlying parameterizations one can visit  
% http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_lifereg_sect019.htm
%
% About the routine:
% First, be informed that John D’Errico’s routines “fminsearchbnd” as well
% as “Adaptive Robust Numerical Differentiation” are required. These can be 
% found and downloaded here:
%
% http://www.mathworks.com/matlabcentral/fileexchange/8277
% http://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation
%
% The first might be used for minimization of the -log-likelihood and 
% the second will be used for the derivation of the hessian matrix. In some special 
% cases the covariance matrix might be derived in closed form, however the
% tools “Adaptive Robust Numerical Differentiation” are accurate and are 
% used throughout the “aft” routine. The folder contains other subroutines 
% but only the “aft” is to be called by the user.
%--------------------------------------------------------------------------
% Input arguments:
%--------------------------------------------------------------------------
% time:        The (possibly right censored) time variable. That is T of the above model.
% Z:           The covariate (design) matrix, each column is one covariate. The column of 
%              ones for the intercept must be included by the user, otherwise no intercept will be considered.
% status:      A vector of length equal to length(time) that
%              takes the values 0 or 1 or -1 for an event a right censored or a left
%              censored observation respectively. That is, one might have
%              available both left and/or right censored data. However, at
%              least some observations must correspond to events.
% distr:       The desired distribution of log(T). Can be EXACTLY set to one of the following: ‘exponential’,
%              ‘weibull’, ‘lognormal’, ‘loglogistic’, ‘gamma’. The last
%              refers to the generalized gamma model. (see also "examples"
%              file for some additional comments).
%
% Note that no missing values (NaNs) are allowed. Please "clean" your data set
% before using "aft".
%--------------------------------------------------------------------------
% Optional Input arguments: (these can either not be reached at all or set as [] in order to continue to the next input argument).
%--------------------------------------------------------------------------
% init:        A row vector of initial values defined by the user. If you do not want 
%              to define this and move on to the next optional input arguments set it as []. 
%              If initial values are not provided by the user then they are derived by OLS.
% minimizer:   Can be set equal to 1 or 2 or 3. That is 1: the fminsearchbnd routine is
%              used for minimization, 2: the fmincon with Largescale set to off and 
%              the ‘spq’ algorithm is employed. 3. the fmincon with Largescale set to 
%              off and the ‘interior-point’ algorithm is employed. If set
%              to [] or not reached at all then minimizer 1 is chosen by default (for no particular reason).
%              If convergence problems occur be sure to explore all three
%              minimizers first. Then you can explore with other initial values
%              in the previous input argument or other options in the next input argument.
% options:     The options that will be used in the optimset of fminsearchbnd or fmincon 
%              depending which minimizer is selected by the previous input argument. That is,
%              if the minimizer is set to 1, then the options will refer to the fminsearchbnd.
%              If the minimizer is set to 2 or 3, then the options will refer to fmincon. The 
%              options set in this argument will replace the default options mentioned in the 
%              minimizer input argument. The default will take place if this input argument is 
%              set to [] or not reached at all. See also MATLAB
%              documentation for the "options" of fmincon or fminsearch. 
%printresults: If set to 1 it prints out a results table in a user friendly format
%              with all the results and the methods employed. If set to any other 
%              value or not reached at all it has no effect. 
%--------------------------------------------------------------------------
% Output arguments:
%--------------------------------------------------------------------------
% pars:        a vector with the estimated parameters. The last element is 
%              the estimated scale parameter. In the case of the "gamma"
%              there is an additional shape parameter. In the case of the
%              "exponential" the scale parameter is held fixed and equal to 1 and 
%              omitted from pars vector.
% covpars:     the variance covariance matrix of the estimated parameters.
% SE:          the standard errors of the estimated regression parameters.
% CI:          95% confidence intervals of the estimated regression parameters.
% Zscores:     Z-scores for the estimated regression parameters
% pvalues:     the p-values for the corresponding regression parameters
% gval:        the achieved maximum of the log-likelihood.
% exitflag:    as defined in the fmincon or fminsearchbnd (depending on the minimizer used)
% hess:        the hessian matrix. This is obtained by “Adaptive Robust
%              Numerical Differentiation” tools provided by John D'errico.
%              If the evaluation of the hessian fails then the one
%              provided by the fmincon is provided and this is noted on top
%              of the results table (if the results table is requested). In
%              that case results might not be accurate and a warning pops
%              up even if a table results is not requested.
% 
% See the “examples” file for a cell by cell illustration and some further
% comments and references.
%
% Generally this code was built as a useful tool, so that interested MATLAB 
% users may fit this kind of models. I do not claim that this code is more efficient
% than other competitive products. It simply builds the corresponding
% -log-likelihoods and uses John D'Erricos fminsearchbnd or built in "fmincon" for minimization.
% The hessian (and the covariance matrices) are (in all cases) calculated 
% numerically using also John D'Erricos “Adaptive Robust Numerical Differentiation” tools.  
% 
% The program was built with version R2010a. Hopefully it would run in
% later (and previous) versions but I am not sure.
% The optimization toolbox is required.
%
%--------------------------------------------------------------------------
% Author:
% Leonidas E. Bantis
% Dept. of Statistics and Actuarial-Financial Mathematics
% University of the Aegean
% Samos Island, Greece.
% e-mail: lbantis@aegean.gr, leobantis@gmail.com
% Release: 1
% Last Updated: September 4th 2012
%--------------------------------------------------------------------------




%----Some error checking:
if size(time,1)==1;time=time';end
if size(status,1)==1;status=status';end
if min(size(status))>1;error('"status" must be a column vector');end
if min(size(time))>1;error('"time" must be a column vector');end
if min(time)<0;error('"time" values must be non-negative');end


if size(Z,1)~=length(time);error('Matrix Z must have the same number of rows as "time"');end

if sum(isnan(time))>0;error('No NaNs are allowed in the data, check "time"');end
if sum(isnan(status))>0;error('No NaNs are allowed in the data, check "status"');end
if sum(sum(isnan(Z)))>0;error('No NaNs are allowed in the data, check "Z"');end

if length(find(status==0))+length(find(status==1))+length(find(status==-1))~=length(time);error('status must be a boolean vector taking values 0 or 1');end



%deal with the correct status coding:
if isempty(find(status==0, 1));error('Analysis cannot be performed if all data are censored');end
%minimizer:
if nargin>=6
    if isempty(minimizer)==1;minimizer=1;end
    if length(minimizer)>1;error('The "minimizer" input argument must be equal to 1 or 2 or 3 or []');end
    if isempty(find(minimizer==[1 2 3], 1));error('The "minimizer" input argument must be equal to 1 or 2 or 3 or []');end
end


if nargin<4;error('Too few input arguments');end

if nargin==4;%distribution names
    distNames = {'exponential', 'weibull', 'loglogistic', 'lognormal', 'gamma'};
    nd = strmatch(lower(distr), distNames);
    if numel(nd) ~= 1
        error('Unknown input argument defined for "distr"');
    else
        distr = distNames{nd};
    end
end

if nargin==4;init=[];minimizer=1;options=[];printresults=0;end
if nargin==5;minimizer=1;options=[];printresults=0;end
if nargin==6;options=[];printresults=0;end
if nargin==7;printresults=0;end
if nargin==8 && printresults~=1;printresults=0;end
if nargin>8;error('Too many input arguments');end

%------End of error checking--
    

switch distr
    
    case 'exponential'
        [pars covpars SE CI Zscores pvalues gval exitflag hess fminconhess]=aftexponential(time,Z,status,init,minimizer,options,printresults);
    case 'weibull'
        %[pars covpars SE CI Zscores pvalues gval exitflag hess fminconhess]=aftweibull(time,Z,status,init,minimizer,options,printresults);
         [pars]=aftweibull2(time,Z,status,init,minimizer,options,printresults);
    case 'gamma'
        [pars covpars SE CI Zscores pvalues gval exitflag hess fminconhess]=aftgamma(time,Z,status,init,minimizer,options,printresults);
    case 'loglogistic'
        [pars covpars SE CI Zscores pvalues gval exitflag hess fminconhess]=aftloglogistic(time,Z,status,init,minimizer,options,printresults);
    case 'lognormal'
        [pars covpars SE CI Zscores pvalues gval exitflag hess fminconhess]=aftlognormal(time,Z,status,init,minimizer,options,printresults);
end

if sum(sum(isnan(hess)))>0 || sum(sum(imag(hess)))>0
warning('warn:hessian','Hessian could not be evaluated at the achieved maximum');
end

if minimizer>1
    if fminconhess==1
        warning('warn:fminconhess','Hessian is obtained by fmincon: results may be inaccurate');
    end
end


end