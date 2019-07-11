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
[X_data, Y_data] = datasetReduction(X, Y, H_est);
X_test = [X(:), Y(:)];

%% Gaussian Process

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = @covSEiso;              % ARD SE
likfunc = @likGauss;              % Gaussian likelihood
% The hyperparameter struct
hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);

tic
% Optimization
hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc,...
    likfunc, X_data, Y_data);
toc

% Let us make a prediction with above parameters!
[mu, s2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc,...
    X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest

%% Approximation result
H_approx = reshape(mu,[m,n]);
err = immse (H, H_approx)

%% PLOT

% % Altitude data
% figure
% mesh(X,Y,H)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% zlim([-50 40])
% saveas(gcf,'origin.png')
% 
% % Noisy data
% figure
% mesh(X,Y,H_est)
% xlabel('x[mm]')
% ylabel('y[mm]')
% zlabel('h[mm]')
% zlim([-50 40])
% saveas(gcf,'noisy.png')

% Approximated data
figure
mesh(X,Y,H_approx)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'approx.png')
