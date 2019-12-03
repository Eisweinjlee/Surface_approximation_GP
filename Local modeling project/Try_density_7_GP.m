%% GP - Try the density of data
% Author: LI, Yang
% Date: Dec 23rd, 2019
close all
clear

%% load the dataset
load local_sparse_7_data.mat
% X_data, Y_data, X_test inside for training
% X_local, Y_local inside for plot

[m,n] = size(X_local); % Grid information

%% Gaussian Process

N = 200;  % iteration times limitation

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;               % inference with Guassian Likelihood
% The hyperparameter struct
hyp_init = struct('mean', [], 'cov', [0 0 0], 'lik', -1);
    
% Optimization hyperparameters  - exact GP
tic;
hyp_exactGP = minimize(hyp_init, @gp, -N, infmethod, meanfunc, covfunc,...
    likfunc, X_data, Y_data);
Time_of_ExactGP = toc;

% The dense prediction
[mu, s2] = gp(hyp_exactGP, infmethod, meanfunc, covfunc, likfunc,...
    X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest

% reshape the result
H_local_exactGP = reshape(mu,[m,n]);

%% Sparse Gaussian Process

% inducing points
% mid = floor(length(X_test)/2);
% xu = X_test(mid-3:mid+3); % (1) bad choice of initial inducing set
xu = X_test(1:25:end,:); % (2) generous choice of initial inducing set
cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s', 0));
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1

tic; hyp_init.xu = xu;
hyp_sparseGP = minimize(hyp_init, @gp, -N, inff, meanfunc, cov, likfunc,...
    X_data, Y_data);
Time_of_SparseGP = toc;

[ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
    X_data, Y_data, X_test);

% reshape the result
H_local_sparse = reshape(ymu,[m,n]);

% % Marginal likelihood and derivatives
% [nlZ,dnlZ] = gp(hyp_sparseGP,infmethod, meanfunc, cov,...
%     likfunc, X_data, Y_data);

%% Plot
Time_of_ExactGP
Time_of_SparseGP
diff_ExactAndSparse = immse(mu,ymu)

figure
mesh(X_local.*170, Y_local.*80, H_local_exactGP)
figure
mesh(X_local.*170, Y_local.*80, H_local_sparse)