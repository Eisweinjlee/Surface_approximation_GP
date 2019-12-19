%% Unit-symmetric Gaussian probability distribution function
% to study about the skewed distribution of Guassian pdf.
% Author: Yang LI
% Date: Dec 16th, 2019
close all
clear

V = 10;
alpha = 0.9;
lambdaX = 10; lambdaY = 10;
x1pos = 0:0.1:5; x1neg = -5:0.1:0; x1neg(end)=[];
x2pos = 0:0.1:5; x2neg = -5:0.1:0; x2neg(end)=[];
x2 = [x2neg,x2pos];

%% 1D situation
Sneg = inv(sqrt(2*pi*lambdaX)) * exp(-0.5*inv(lambdaX)*(alpha/(2*alpha-1)*x1neg).^2);
Spos = inv(sqrt(2*pi*lambdaX)) * exp(-0.5*inv(lambdaX)*(alpha*x1pos).^2);
S = [Sneg,Spos];
figure;plot([x1neg,x1pos],S)

%% 2D situation
[Xpos,Ypos] = meshgrid(x1pos,x2);
[Xneg,Yneg] = meshgrid(x1neg,x2);
Sigma = [lambdaX,0;0,lambdaY];

Hpos = V/(2*pi*det(Sigma)^0.5)*exp(-0.5*(inv(lambdaX)*(alpha*Xpos).^2 + ...
    inv(lambdaY)*Ypos.^2));
Hneg = V/(2*pi*det(Sigma)^0.5)*exp(-0.5*(inv(lambdaX)*(alpha/(2*alpha-1)*Xneg).^2 + ...
    inv(lambdaY)*Yneg.^2));

% the wall effect
xf = 5; yr = -5; yl = 5;
KK = V/(2*pi*det(Sigma)^0.5);

Hpos = Hpos + KK *exp(-0.5*(inv(lambdaX)*(alpha*Xpos).^2 + ...
    inv(lambdaY)*(-(Ypos-(2*yl-0))).^2)); % yl - pos
Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha/(2*alpha-1)*Xneg).^2 + ...
    inv(lambdaY)*(-(Yneg-(2*yl-0))).^2)); % yl - neg

Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*Xpos).^2 + ...
    inv(lambdaY)*(-(Ypos-(2*yr-0))).^2)); % yr - pos
Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha/(2*alpha-1)*Xneg).^2 + ...
    inv(lambdaY)*(-(Yneg-(2*yr-0))).^2)); % yr - neg

Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xneg-(2*xf-0)))).^2 + ...
    inv(lambdaY)*Yneg.^2)); % xf - neg
Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xneg-(2*xf-0)))).^2 + ...
    inv(lambdaY)*(-(Yneg-(2*yl-0))).^2)); % xf & yl - neg
Hneg = Hneg + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xneg-(2*xf-0)))).^2 + ...
    inv(lambdaY)*(-(Yneg-(2*yr-0))).^2)); % xf & yr - neg

Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xpos-(2*xf-0)))).^2 + ...
    inv(lambdaY)*Ypos.^2)); % xf - pos
Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xpos-(2*xf-0)))).^2 + ...
    inv(lambdaY)*(-(Ypos-(2*yl-0))).^2)); % xf & yl - pos
Hpos = Hpos + KK*exp(-0.5*(inv(lambdaX)*(alpha*(-(Xpos-(2*xf-0)))).^2 + ...
    inv(lambdaY)*(-(Ypos-(2*yr-0))).^2)); % xf & yr - pos

H = [Hneg,Hpos];

X = [Xneg,Xpos]; Y = [Yneg,Ypos];
figure;mesh(X,Y,H)