%% Comparison on FITC and VFE demo
% LI Yang// Oct 22nd, 2019
clear
close all

% Method selection: 0: VFE, 1: FTIC, 2: full GP
selection_flag = 0;

%% The data generation

% latent function
x = linspace(-2,2);
y = 0.6*sin(2*pi.*x + 0.2) + 0.8*cos (pi.*x -0.5);
% plot(x,y)
% axis([-2 2 -2 2])

[m,n] = size(y);

% Gaussian noise
var = 0.2;
y1 = y + var .* randn(m,n);
y2 = y + var .* randn(m,n);
y3 = y + var .* randn(m,n);
y4 = y + var .* randn(m,n);

% Cauchy noise
cau = 0.08;
y5 = y + cau * tan(pi*(rand(m,n)-1/2));

%% Dataset

X_data = [x';x';x';x';x'];
Y_data = [y1';y2';y3';y4';y5'];

% plot(X_data,Y_data,'*')
% axis([-2 2 -2 2])

X_test = x';

%% Gaussian Process

N = 200;  % iteration times limitation

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;               % inference with Guassian Likelihood
% The hyperparameter struct
hyp_init = struct('mean', [], 'cov', [0 0], 'lik', -1);

if selection_flag == 2 % full GP regression
    
    % Optimization hyperparameters
    tic
    hyp_full = minimize(hyp_init, @gp, -N, infmethod, meanfunc, covfunc,...
        likfunc, X_data, Y_data);
    toc
    
    % The dense prediction
    [mu, s2] = gp(hyp_full, infmethod, meanfunc, covfunc, likfunc,...
        X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest 
    
elseif selection_flag == 0 || selection_flag == 1 % Sparse approximation
    
    % inducing points
    xu = X_test(45:55); % (1) bad choice of initial inducing set
    % xu = X_test(1:2:end); % (2) generous choice of initial inducing set
    cov = {'apxSparse', covfunc, xu};
    inff = @(varargin) infmethod(varargin{:},struct('s', selection_flag));
    % VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1
    
    tic
    hyp_init.xu = xu;
    hyp = minimize(hyp_init, @gp, -N, inff, meanfunc, cov, likfunc,...
        X_data, Y_data);
    toc
    
    [ymu,ys2] = gp(hyp, inff, meanfunc, cov, likfunc,...
        X_data, Y_data, X_test);
    
    % Marginal likelihood and derivatives
    [nlZ,dnlZ] = gp(hyp,infmethod, meanfunc, cov,...
        likfunc, X_data, Y_data);
    
end


%% Plot
if selection_flag == 2
    
    % Analysis
    error_mse = immse(mu,y')
    
    f = [mu+2*sqrt(s2); flipdim(mu-2*sqrt(s2),1)];
    
    figure
    hold on
    fill([X_test; flipdim(X_test,1)], f, [0.9 0.9 0.9]);
    
    plot(X_data,Y_data,'b*')
    plot(x,y,'--b','Linewidth',1.5)
    plot(X_test,mu,'r','Linewidth',1.5)
    
    hold off
    axis([-2 2 -2 2])
    grid on
    
elseif selection_flag == 0 || selection_flag == 1
    
    % Analysis
    error_mse = immse(ymu,y')
    
    f = [ymu+2*sqrt(ys2); flipdim(ymu-2*sqrt(ys2),1)];
    
    figure
    hold on
    fill([X_test; flipdim(X_test,1)], f, [0.9 0.9 0.9]);
    
    plot(X_data,Y_data,'b*')
    plot(x,y,'--b','Linewidth',1.5)
    plot(X_test,ymu,'r','Linewidth',1.5)
    
    plot(xu,2*ones(length(xu),1),'o')
    plot(hyp.xu,-2*ones(length(hyp.xu),1),'o')
    hold off
    axis([-2 2 -2 2])
    grid on
    
end