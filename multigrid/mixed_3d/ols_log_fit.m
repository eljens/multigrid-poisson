function [alpha,beta] = ols_log_fit(err,h)
    err = err(:);
    h = h(:);
    X = [ones(size(h)),log(h)];
    parm = (X'*X) \ (X'*log(err));
    alpha = parm(1);
    beta = parm(2);
end
