%% |mxe_getsolverhandle|
% *Note:* This is a private function.
%
% Returns function handle to a manopt solver given its name
%
% *Syntax*
%
%   h = mxe_getsolverhandle(solverName)
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function h = mxe_getsolverhandle(solverName)

    switch lower(solverName)
        case {'cg', 'conjugategradient'}
            h = @conjugategradient;
        case 'lbfgs'
            h = @lbfgs;
        case {'tr', 'trustregions'}
            h = @trustregions;
        case {'sd', 'gd', 'steepestdescent', 'gradientdescent'}
            h = @steepestdescent; 
        case {'sgd', 'stochasticgradientdescent'}
            h = @sgd; 
        case 'lbfgs_transportfree'
            h = @lbfgs_TransportFree;
        case 'lbfgs_transportfreecholesky'
            h = @lbfgs_TransportFreeCholesky;
        case 'lbfgs_mixest'
            h = @lbfgs_MIXEST;
        case 'lbfgs_manopt'
            h = @lbfgs_MANOPT;
        %TODO
        otherwise
            error('Solver "%s" not defined in mxe_getsolverhandle.', solverName)
    end
end

