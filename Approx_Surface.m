clear
close all

%% Vessel parameter initailize
load('Vessel_XY.mat') % The meshgrid: X & Y

%% 2D altitude data with 94 x 100 elements
load('surfacedata.mat')
H = surfacedata;
clear surfacedata
m = length(H(:,1));
n = length(H(1,:));

%% The Dataset Generation
[X_latent, Y_latent] = datasetReduction(X, Y, H);

s = 5;  % variance for Observation
s_in = 5e-1;  % variance for Input

X_data = [];
Y_data = [];
% Gaussian Noise
for i = 1:1:10
X_noise = X_latent + s_in * randn(length(X_latent(:,1)),length(X_latent(1,:)));
Y_noise = Y_latent + s * randn(length(Y_latent(:,1)),length(Y_latent(1,:)));

X_data = [X_data;X_noise];
Y_data = [Y_data;Y_noise];

end

X_test = [X(:), Y(:)];

%% Gaussian Process

N = 100;  % iteration times limitation

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;               % inference with Guassian Likelihood
% The hyperparameter struct
hyp_init = struct('mean', [], 'cov', [0 0 0], 'lik', -1);

% % Optimization hyperparameters
% tic
% hyp_full = minimize(hyp_init, @gp, -N, infmethod, meanfunc, covfunc,...
%     likfunc, X_data, Y_data);
% toc
% 
% % The dense prediction
% [mu, s2] = gp(hyp_full, infmethod, meanfunc, covfunc, likfunc,...
%     X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest

%% Sparse approximation

% inducing points
xu = X_test(1:15:end,:); cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s',0.0));  
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1

tic
hyp_init.xu = xu;
hyp = minimize(hyp_init, @gp, -75, inff, meanfunc, cov, likfunc,...
    X_data, Y_data);
toc

[ymu,ys2] = gp(hyp, inff, meanfunc, cov, likfunc,...
    X_data, Y_data, X_test);

% % Marginal likelihood and derivatives
% [nlZ,dnlZ] = gp(hyp,infmethod, meanfunc, cov,...
%     likfunc, X_data, Y_data)


%% Approximation result
% % H_approx = reshape(mu,[m,n]);
% err1 = immse (H, H_approx)

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


% % Approximated data
% figure
% mesh(X,Y,H_approx)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% zlim([-50 40])
% % saveas(gcf,'approx.png')

% Sparse approximated
figure
mesh(X,Y,H_sparse)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'sparse.png')

% inducing point optimization
figure
plot(xu(:,1),xu(:,2),'bx','LineWidth',0.8)
hold on
plot(hyp.xu(:,1),hyp.xu(:,2),'rx','LineWidth',0.8)
