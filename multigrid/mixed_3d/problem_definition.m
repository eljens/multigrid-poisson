function [ufun,ffun,dxfun,dyfun] = problem_definition()
    % Function periods
    kx = 4;
    ky = 5;
    kz = 6;

    % Test function
    ufun = @(x,y,z)(sin(kx*x).*sin(ky*y).*sin(kz*z));
    
    % x derivative
    dxfun = @(x,y,z)(kx*cos(kx*x).*sin(ky*y).*sin(kz*z));
    
    % y derivative;
    dyfun = @(x,y,z)(ky*sin(kx*x).*cos(ky*y).*sin(kz*z));
    
    % Laplacian of the test function
    ffun = @(x,y,z)(-(kx^2+ky^2+kz^2)*sin(kx*x).*sin(ky*y).*sin(kz*z));
end