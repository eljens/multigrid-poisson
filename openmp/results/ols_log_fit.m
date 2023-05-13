function [alpha,beta] = ols_log_fit(err,h)
%OLS_LOG_FIT
% This function computes the slope beta and the intercept alpha of a
% straight line trough the errors in log space using ordinary least squares. 
% This function can be used to verify if the order of convergence 
% implementation matches the theoretical order.
%
% Syntax: [alpha,beta] = ols_log_fit(err,h)
%
% Inputs:
%   err: The errors measured in some norm for spacing h. 1D or 2D matrix.
%   h:   The equidistant grid spacings. 1D or 2D matrix.
%
% Outputs:
%   alpha: Intercept. Scalar.
%   beta:  Slope. Scalar
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 18th of November 2022
    % Making sure that the dimension of the input is the same
    err = err(:);
    h = h(:);
    assert(size(h,1) == size(err,1),"The inputs must have the same length");
    X = [ones(size(h)),log(h)];
    % Solving for the parameters
    parm = (X'*X) \ (X'*log(err));
    % Extracting the coefficients
    alpha = parm(1);
    beta = parm(2);
end

