function [ufun,ffun,dxfun,dyfun] = problem_definition()
%PROBLEM_DEFINITION
% This function returns four functions computing the exact solution, the
% right hand side, and the Neumann conditions respectively.
%
% Syntax: [ufun,ffun,dxfun,dyfun] = problem_definition()
%
% Outputs:
%   ufun: The true solution/Dirichlet condition.
%   ffun: The right hand side.
%   dxfun: The x derivative of u(x,y,z).
%   dyfun: The y derivative of u(x,y,z).
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 18th of November 2022
    % Function periods
    kx = 5.56;
    ky = 19;
    kz = -0.34;

    % Test function
    ufun = @(x,y,z)(sin(kx*x).*sin(ky*y).*sin(kz*z));
    
    % x derivative
    dxfun = @(x,y,z)(kx*cos(kx*x).*sin(ky*y).*sin(kz*z));
    
    % y derivative;
    dyfun = @(x,y,z)(ky*sin(kx*x).*cos(ky*y).*sin(kz*z));
    
    % Laplacian of the test function
    ffun = @(x,y,z)(-(kx^2+ky^2+kz^2)*sin(kx*x).*sin(ky*y).*sin(kz*z));
end