clear
close all

%% Vessel parameter initailize
load('Vessel_XY.mat') % The meshgrid

%% 2D altitude data with 94 x 100 elements
load('surfacedata.mat')
H = surfacedata;
clear surfacedata

%% Noisy data
s = 5;  % variance

m = length(H(:,1));
n = length(H(1,:));
H_est = H + s * randn(m,n);   % Gaussian noise
clear m n

%% Fit with polynomial model
x = X(:);
y = Y(:);
h = H_est(:);

f = fit([x,y],h,'poly55');
% plot(f,[x,y],h)

%% The approximate model
H_approx = f(X,Y);
err = immse(H, H_approx)

%% PLOT

% Altitude data
figure
mesh(X,Y,H)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'origin.png')

% Noisy data
figure
mesh(X,Y,H_est)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'noisy.png')

% Approximated data
figure
mesh(X,Y,H_approx)
xlabel('x[mm]')
ylabel('y[mm]')
zlabel('h[mm]')
zlim([-50 40])
% saveas(gcf,'approx.png')