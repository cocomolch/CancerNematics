%%fractional derivative Kelvin-Voigt-model
%%(special case of fractional KV model)
%%Philip Friedrich 23.03.2023

    %%%     https://doi.org/10.1039/D0SM00354A
    %%%     https://doi.org/10.1371/journal.pone.0143090
    %%%     10.1016/j.scitotenv.2017.09.206

function [mu,alpha,mu_z,mu_imag,resnorm,residual] = frac_derivative_KV_model(f,G1,G2)
    % define the model as an 'anonymous' function
    model = @(x,xdata)(x(1)^(1-x(2)))*(1i*2*pi*xdata).^(x(2)) + x(3) + 1i.*x(4);
    
    % check if G', G'', or both are provided
    if exist('G1','var') && exist('G2','var')
        fun = @(x,xdata) [real(model(x,xdata)); imag(model(x,xdata))];
    elseif exist('G1','var') && ~exist('G2','var')
        fun = @(x,xdata) real(model(x,xdata));
    elseif ~exist('G1','var') && exist('G2','var')
        fun = @(x,xdata) imag(model(x,xdata));
    end

    % define starting values for fit parameters  
    x_0 = [1000, 0.3, 1000, 100];
    
    % define lower and upper boundary
    lb=[0, 0, 0, 0];
    ub=[1000000, 1, 1000000, 100000];

    % get frequency data
    xdata = f(~isnan(f));

    % get G' and G'' data
    if exist('G1','var') && exist('G2','var')
        ydata = [G1(~isnan(G1)); G2(~isnan(G2))];
    elseif exist('G1','var') && ~exist('G2','var')
        ydata = G1(~isnan(G1));
    elseif ~exist('G1','var') && exist('G2','var')
        ydata = G2(~isnan(G2));
    end

    % suppress lsqcurvefit output
    options = optimset('Display','off');

    % fit function
    [x,resnorm,residual] = lsqcurvefit(fun,x_0,xdata,ydata,lb,ub,options);

    % define output variables
    mu = x(1);
    alpha = x(2);
    mu_z = x(3);
    mu_imag = x(4);
end
