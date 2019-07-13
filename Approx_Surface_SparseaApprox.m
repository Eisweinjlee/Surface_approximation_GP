clear
close all

%% Vessel parameter initailize
load('Vessel_XY.mat') % The meshgrid

%% 2D altitude data with 94 x 100 elements
load('surfacedata.mat')
H = surfacedata;
clear surfacedata
m = length(H(:,1));
n = length(H(1,:));

%% Noisy data
s = 5;  % variance
H_est = H + s * randn(m,n);   % Gaussian noise

%% The dataset
% simply reduce to 1/4 data, by jump 2 steps

% Naive reduction
[X_data, Y_data] = datasetReduction(X, Y, H_est);

% % full dataset
% X_data = [X(:), Y(:)];
% Y_data = H_est(:);

X_test = [X(:), Y(:)];

%% Gaussian Process

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;               % inference with Guassian Likelihood
% The hyperparameter struct
hyp_init = struct('mean', [], 'cov', [0 0 0], 'lik', -1);

% Optimization hyperparameters
tic
hyp_full = minimize(hyp_init, @gp, -200, infmethod, meanfunc, covfunc,...
    likfunc, X_data, Y_data);
toc

% The dense prediction
[mu, s2] = gp(hyp_full, infmethod, meanfunc, covfunc, likfunc,...
    X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest

% % Marginal likelihood and derivatives
% [nlZ,dnlZ] = gp(hyp_full,infmethod, meanfunc, covfunc,...
%     likfunc, X_data, Y_data);

%% Sparse approximation

% inducing points
xu = X_test(1:10:end,:); cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s',1.0));  
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1

tic
[ymu,ys2] = gp(hyp_full, inff, meanfunc, cov, likfunc,...
    X_data, Y_data, X_test);
toc

%% Approximation result
H_approx = reshape(mu,[m,n]);
err1 = immse (H, H_approx)

H_sparse = reshape(ymu,[m,n]);
err2 = immse (H, H_sparse)


%% PLOT

% % Altitude data
% figure
% mesh(X,Y,H)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% zlim([-50 40])
% % saveas(gcf,'origin.png')

% % Noisy data
% figure
% mesh(X,Y,H_est)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% zlim([-50 40])
% % saveas(gcf,'noisy.png')

% Approximated data
figure
mesh(X,Y,H_approx)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'approx.png')

% Sparse approximated
figure
mesh(X,Y,H_sparse)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'sparse.png')