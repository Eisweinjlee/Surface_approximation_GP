% Find the optimal Sigma and Volumn factor
% Original version date: Oct 8th, 2019
% Date: Dec 5th, 2019
% Author: Yang LI
close all
clear

%% 1. Initialization
excavator_data;
docName = "Approx_Surface\ErrorData training project\dataset_20191108_normalized\";

% load the empty vessel
load(docName + "00.mat",'dep')
H0 = dep;
[m,n] = size(dep);
H_data = zeros(m,n,61);

% load others
number = 0;
for i = 1:9
    for j = 1:5
        number = number + 1;
        filename = docName + num2str(10*i+j)+".mat";
        load(filename)
        H_data(:,:,number) = dep - H0;
    end
end
for i = 1:16
    number = number + 1;
    filename = docName + num2str(100+i)+".mat";
    load(filename)
    H_data(:,:,number) = dep - H0;
end

%% 2. The loading center data
Xc = [25*ones(5,1);50*ones(5,1);75*ones(15,1);50*ones(10,1);25*ones(10,1)];
Yc = [zeros(15,1);40*ones(5,1);-40*ones(10,1);40*ones(10,1);-40*ones(5,1)];
Xc = [Xc;62.5;37.5;25;25;37.5;62.5;75;75;50;50;62.5;37.5;37.5;37.5;62.5;62.5]*1.7;
Yc = [Yc;40;40;20;-20;-40;-40;-20;20;20;-20;0;0;-20;20;20;-20];

parameters = zeros(3,61);
for i = 1:number
%% optimization program initialization
lambdaX = 1709.10978; lambdaY = 2274.09987; kV = 3.8;

depH = H_data(:,:,i);    
Sigma = [lambdaX,0; 0,lambdaY];
depx = Xc(i); depy = Yc(i); c = [depx, depy]; 
the = atan((depy-Pe(2))/(depx-Pe(1)));
V = sum(depH ,'all'); % 1.8e+5 is the actual volume

ModelH = function_input_2d(X,Y,c,kV*V,Sigma,the,xf,yr,yl);

% figure; subplot(1,2,1); mesh(X,Y,ModelH);xlim([0 170]);zlim([0 40]);
% subplot(1,2,2); mesh(X,Y,depH);xlim([0 170]);zlim([0 40]);

%% Error analysis
error_before = immse(depH, ModelH);

%% Optimization
theta0 = [2000,2000,3];
lb = [1000,1000,2.5];
ub = [3000,3000,4];

fun = @(theta)immse(depH, function_input_2d(X,Y,c,theta(3)*V,[theta(1),0;0,theta(2)],the,xf,yr,yl));

options = optimset('Display','iter','PlotFcns',@optimplotfval);

tic
theta = fmincon(fun,theta0,[],[],[],[],lb,ub);
% theta = fminsearch(fun,theta0);
toc

parameters(:,i) = theta';

%% result
lambdaX = theta(1); lambdaY = theta(2); kV = theta(3);

Sigma = [lambdaX,0; 0,lambdaY];
ModelH = function_input_2d(X,Y,c,kV*V,Sigma,the,xf,yr,yl);

% error_before
error_after = immse(depH, ModelH);
error_before - error_after

figure
subplot(1,2,1)
mesh(X,Y,ModelH)
xlim([0 170])
zlim([0 40])

subplot(1,2,2)
mesh(X,Y,depH)
xlim([0 170])
zlim([0 40])
% 
% input('Next data?[Enter]')
% close all
end

%% summary
maybeGoodResult = mean(parameters,2)
% % resultant parameters are obtained from 7 data

lambdaX = maybeGoodResult(1); lambdaY = maybeGoodResult(2);
kV = maybeGoodResult(3);
% save nominal_modle_para lambdaX lambdaY kV