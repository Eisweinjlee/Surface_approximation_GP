clear
close all

%% Vessel parameter initailize
load('Vessel_XY.mat') % The meshgrid

%% 2D altitude data with 94 x 100 elements
load('surfacedata.mat')
H = surfacedata;
clear surfacedata

%% Noisy data
s = 5e-2;  % variance

m = length(H(:,1));
n = length(H(1,:));
H_est = H + s * randn(m,n);   % Gaussian noise

%% The dataset
x = X(:);
y = Y(:);
h = H_est(:);

X_data = [x,y];
Y_data = h;

X_test = [x,y];

%% Gaussian Process

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = @covSEard;              % ARD SE
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
    X_data, Y_data, X_data);  % Xtrain,Ytrain,Xtest

%% Approximation result
H_approx = reshape(mu,[m,n]);
err = immse (H, H_approx)

%% PLOT

% Altitude data
figure
mesh(X,Y,H)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
saveas(gcf,'origin.png')

% Noisy data
figure
mesh(X,Y,H_est)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
saveas(gcf,'noisy.png')

% Approximated data
figure
mesh(X,Y,H_approx)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
saveas(gcf,'approx.png')
