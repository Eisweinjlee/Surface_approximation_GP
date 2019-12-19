%% Unit-symmetric Gpdf soil loading model
% to use the unit-symmetric Gaussian pdf as soil loading model.
% Author: Yang LI
% Date: Dec 19th, 2019

% % [test] input
% clear
% close all
% 
% load("C:\Users\LI\Desktop\Approx_Surface\Vessel_XY.mat")
% c = [85,0];
% V = 55000;
% Sigma = [2000,0;0,2000];
% alpha = 0.85;

function H = unit_sym_input_2d(X,Y,c,V,Sigma,alpha)

% Inputs:
% X,Y - the 2D mesh grid, size: two 94x100 matrices
% c - the loading center, size: 2d vector
% V - the loading volume, size: real number
% Sigma - the covariance matrix, size: 2x2 matrix
% alpha - the parameter of skewness, size: real number

xc = c(1); yc = c(2);
[~,I] = min(abs(X(1,:)-xc));

Xpos = X(:,I+1:end); Xneg = X(:,1:I);
Ypos = Y(:,I+1:end); Yneg = Y(:,1:I);

lambdaX = Sigma(1,1); lambdaY = Sigma(2,2);
KK = V/(2*pi*det(Sigma)^0.5);

Hpos = KK*exp(-0.5*(inv(lambdaX)*(alpha*(Xpos-xc)).^2 + ...
    inv(lambdaY)*(Ypos-yc).^2));
Hneg = KK*exp(-0.5*(inv(lambdaX)*(alpha/(2*alpha-1)*(Xneg-xc)).^2 + ...
    inv(lambdaY)*(Yneg-yc).^2));

% % wall effect
% xf = 170; yr = -80; yl = 80; % the position of walls
% 
% Hpos = Hpos + KK *exp(-0.5*(inv(lambdaX)*(alpha*(Xpos-xc)).^2 + ...
%     inv(lambdaY)*(-(Ypos-(2*yl-yc))).^2)); % yl - pos
% Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha/(2*alpha-1)*(Xneg-xc)).^2 + ...
%     inv(lambdaY)*(-(Yneg-(2*yl-yc))).^2)); % yl - neg
% 
% Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(Xpos-xc)).^2 + ...
%     inv(lambdaY)*(-(Ypos-(2*yr-yc))).^2)); % yr - pos
% Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha/(2*alpha-1)*(Xneg-xc)).^2 + ...
%     inv(lambdaY)*(-(Yneg-(2*yr-yc))).^2)); % yr - neg
% 
% Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xneg-(2*xf-xc)))).^2 + ...
%     inv(lambdaY)*(Yneg-yc).^2)); % xf - neg
% Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xneg-(2*xf-xc)))).^2 + ...
%     inv(lambdaY)*(-(Yneg-(2*yl-yc))).^2)); % xf & yl - neg
% Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xneg-(2*xf-xc)))).^2 + ...
%     inv(lambdaY)*(-(Yneg-(2*yr-yc))).^2)); % xf & yr - neg
% 
% Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xpos-(2*xf-xc)))).^2 + ...
%     inv(lambdaY)*(Ypos-yc).^2)); % xf - pos
% Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xpos-(2*xf-xc)))).^2 + ...
%     inv(lambdaY)*(-(Ypos-(2*yl-yc))).^2)); % xf & yl - pos
% Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xpos-(2*xf-xc)))).^2 + ...
%     inv(lambdaY)*(-(Ypos-(2*yr-yc))).^2)); % xf & yr - pos

H = [Hneg,Hpos];
end

% % [test] plot
% X = [Xneg,Xpos]; Y = [Yneg,Ypos];
% figure;mesh(X,Y,H); xlim([0 170]);