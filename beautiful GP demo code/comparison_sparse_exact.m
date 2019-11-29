%% Comparison on Sparse and exact GP
% LI Yang// Nov 7th, 2019
clear
close all

% Method selection: (0)VFE, (1)FTIC
selection_flag = 0;

% New data(0) or load data(1)
data_type_flag = 1;

if data_type_flag == 0
%% The data generation
    
    % latent function
    x = linspace(-2,2,1000);
    y = 0.6*sin(2*pi.*x + 0.2) + 0.8*cos (pi.*x -0.5);
    % plot(x,y)
    % axis([-2 2 -2 2])
    
    [m,n] = size(y);
    
    % Gaussian noise
    var = 0.25;
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
    
%     save('comparison_sparse_exact_data.mat','X_data','Y_data','X_test','x','y');
%     save('comparison_very_dense_data.mat','X_data','Y_data','X_test','x','y');
elseif data_type_flag == 1
% load data
    load('comparison_sparse_exact_data.mat')
%     load('comparison_very_dense_data.mat')
end

%% Gaussian Process

N = 200;  % iteration times limitation

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;               % inference with Guassian Likelihood
% The hyperparameter struct
hyp_init = struct('mean', [], 'cov', [0 0], 'lik', -1);
    
% Optimization hyperparameters
tic;
hyp_exactGP = minimize(hyp_init, @gp, -N, infmethod, meanfunc, covfunc,...
    likfunc, X_data, Y_data);
Time_of_ExactGP = toc;

% The dense prediction
[mu, s2] = gp(hyp_exactGP, infmethod, meanfunc, covfunc, likfunc,...
    X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest


% inducing points
% mid = floor(length(X_test)/2);
% xu = X_test(mid-3:mid+3); % (1) bad choice of initial inducing set
xu = X_test(1:30:end); % (2) generous choice of initial inducing set
cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s', selection_flag));
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1

tic; hyp_init.xu = xu;
hyp_sparseGP = minimize(hyp_init, @gp, -N, inff, meanfunc, cov, likfunc,...
    X_data, Y_data);
Time_of_SparseGP = toc;

[ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
    X_data, Y_data, X_test);

% % Marginal likelihood and derivatives
% [nlZ,dnlZ] = gp(hyp_sparseGP,infmethod, meanfunc, cov,...
%     likfunc, X_data, Y_data);



%% Plot
    
% Analysis
Time_of_ExactGP
Time_of_SparseGP

exactGP_mse = immse(mu,y')
f1 = [mu+2*sqrt(s2), flipdim(mu-2*sqrt(s2),1)];

sparseGP_mse = immse(ymu,y')
f2 = [ymu+2*sqrt(ys2), flipdim(ymu-2*sqrt(ys2),1)];

figure
hold on
% Dataset
plot(X_data,Y_data,'LineStyle','none','Marker','+','Color','[0.8500, 0.3250, 0.0980]')
plot(x,y,'Linewidth',1,'Color','#77AC30')

% Predictive mean plot
plot(X_test,mu,'--b','Linewidth',1.5,'Color','[0.6350, 0.0780, 0.1840]')
plot(X_test,ymu,'r','Linewidth',1.5,'Color','[0, 0.4470, 0.7410]')

% exact GP predictive variance
% fill([X_test; flipdim(X_test,1)], f1(:), [0.9 0.9 0.9]);
plot(X_test,f1(:,1),'--','Linewidth',1,'Color','[0.6350, 0.0780, 0.1840]')
plot(flipdim(X_test,1),f1(:,2),'--','Linewidth',1,'Color','[0.6350, 0.0780, 0.1840]')

% sparse GP predictive variance
% fill([X_test; flipdim(X_test,1)], f2(:), [0.9 0.9 0.9]);
plot(X_test,f2(:,1),'--','Linewidth',1,'Color','[0, 0.4470, 0.7410]')
plot(flipdim(X_test,1),f2(:,2),'--','Linewidth',1,'Color','[0, 0.4470, 0.7410]')

% Inducing point displacement
plot(xu,2.5*ones(length(xu),1),'r^') % initial
plot(hyp_sparseGP.xu,-2.5*ones(length(hyp_sparseGP.xu),1),'rv') % result

hold off
grid on
axis([-2 2 -2.5 2.5])

titlename = "Number of data: " + size(X_data,1) + ", inducing points: " +...
    size(xu,1);
title(titlename);

% figure
% K = feval(covfunc{:}, hyp_exactGP.cov, sort(X_data));
% contour(K,'Fill','on')