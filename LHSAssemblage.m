function [is_LHS_assemblage,siglam] = LHSAssemblage(sigma,varargin)
%LHSAssemblage Determines whether an assemblage has an LHS model or not
%  This function has one required argument:
%   sigma: a 4-D array, containing the members of the assemblage. The first 
%   two dimensions contain the (unnormalised) quantum states, while the
%   remaining two dimensions are (a,x), such that sigma(:,:,a,x) =
%   \sigma_a|x. 
%
%  is_LHS_assemblage = LHSAssemblage(sigma) is the indicator function for
%  LHS assemblages. It returns 1 if the assemblage is a valid unnormalised
%  LHS assemblage, and 0 otherwise.
%
%  [is_LHS_assemblage,siglam] = LHSAssemblage(sigma) also returns the LHS
%  model which reproduces the assemblage sigma. In particular, siglam is
%  a dB x dB x ndet array, containing the ndet = oa^ma members of the 
%  LHS model, one corresponding to each deterministic strategy for Alice.
%
% This function has one optional argument:
%   nm: (default 0)
%
% is_LHS_assemblage = LHSAssemblage(sigma,nm) is the indicator function for
% LHS assemblages that are additionally required to be normalised when nm =
% 1. In this case, trace(sum(sigma(:,:,:,x),3)) = 1.
%
% is_LHS_assemblage = LHSAssemblage(sigma,nm) can also be used inside CVX
% as a partially specified problem, to enforce the constraint that the CVX
% variable sigma should be an unnormalised LHS assemblage if nm = 0, and a
% normalised LHS assemblage otherwise.
%
% EXAMPLE:
%   cvx_begin
%
%       variable sigma(dB,dB,oa,ma)
%
%       subject to
%
%       LHSAssemblage(sigma) == 1
%
%   cvx_end
%
%  Inside CVX, sigma is an assemblage with ma inputs and oa outcomes for
%  Alice, preparing quantum states of dimension dB x dB for Bob. CVX
%  enforces that sigma should have an LHS model.
%
%   requires: CVX (http://cvxr.com/cvx/), QETLAB (http://www.qetlab.com)
%   authors: Paul Skrzypczyk, Daniel Cavalcanti
%   last updated: March 17, 2016

[nm] = opt_args({0},varargin{:}); 
%if unspecified, it is assumed that the assemblage is unnormalised.

[dB,~,oa,ma] = size(sigma);
% dB = dim. of Bob, oa = # outcomes for Alice, ma = # inputs for Alice

Ndet = oa^ma; % number of deterministic strategies for Alice
SingleParty = genSinglePartyArray(oa,ma); % generate array containing 
                                            %single party strategies 

if isa(sigma,'cvx') == 0 % if sigma isn't a CVX variable  
    % check that the assemblage is valid
    if  NSAssemblage(sigma,nm) == 0
        error('assemblage is not valid')
    end 
end

cvx_begin sdp quiet
    
    variable siglam(dB,dB,Ndet) hermitian semidefinite
    % siglam are the members of the LHS model

    subject to
    
    sigma == squeeze(sum(repmat(siglam,[1,1,1,oa,ma])...
        .*permute(repmat(SingleParty,[1,1,1,dB,dB]),[4,5,3,1,2]),3));
    % sig_a|x == \sum_lam D(a|x,lam) sig_lam 
    
    if nm == 1 % if required to be normalised, sum_lam tr(sig_lam) == 1
        trace(sum(siglam,3)) == 1;
    end

cvx_end

% CVX will return +inf if the problem is infeasible, and 0 if feasible
% this maps {+inf,0} to {0,1}
is_LHS_assemblage = 1-min(cvx_optval,1);

% if there is no LHS model, then return siglam as the empty array
if is_LHS_assemblage == 0
    siglam = [];
end
    
end


    