%% Find the optimal parameters for Gaussian PDF
% Sigma, center move, volume factor and alpha for skewing
% Original version date: Oct 8th, 2019
% Date: Dec 19th, 2019
% Author: Yang LI
close all
clear

%% 1. Initialization
excavator_data;

docName1 = "Approx_Surface\ErrorData training project\dataset_20191108_normalized\";
docName2 = "Approx_Surface\ErrorData training project\dataset_20191209_normalized\";
m = 94; n = 100;
H_data = zeros(m,n,61+36);

% load the empty vessel
load(docName1 + "00.mat",'dep')
H0 = dep;

% load docName1 - Position difference
number = 0;
for i = 1:9
    for j = 1:5
        number = number + 1;
        filename = docName1 + num2str(10*i+j)+".mat";
        load(filename)
        H_data(:,:,number) = dep - H0;
    end
end
for i = 1:16
    number = number + 1;
    filename = docName1 + num2str(100+i)+".mat";
    load(filename)
    H_data(:,:,number) = dep - H0;
end
% load docName2 - Volume difference
for i = 1:9
    for j = 1:4
        number = number + 1;
        filename = docName2 + num2str(10*i+j)+".mat";
        load(filename)
        H_data(:,:,number) = dep - H0;
    end
end

%% 2. The loading center data
% docName1
Xc = [25*ones(5,1);50*ones(5,1);75*ones(15,1);50*ones(10,1);25*ones(10,1)];
Yc = [zeros(15,1);40*ones(5,1);-40*ones(10,1);40*ones(10,1);-40*ones(5,1)];
Xc = [Xc;62.5;37.5;25;25;37.5;62.5;75;75;50;50;62.5;37.5;37.5;37.5;62.5;62.5];
Yc = [Yc;40;40;20;-20;-40;-40;-20;20;20;-20;0;0;-20;20;20;-20];
% docName2
Xc = [Xc;25*ones(4,1);50*ones(4,1);75*ones(12,1);50*ones(8,1);25*ones(8,1)]*1.7;
Yc = [Yc;zeros(12,1);40*ones(4,1);-40*ones(8,1);40*ones(8,1);-40*ones(4,1)];

%% 3. Volume data
Vol_data = ones(61,1);
for i = 1:9
    Vol_data = [Vol_data;0.25;0.50;0.50;0.75];
end

%% 4. optimization program initialization
parameters = zeros(6,number);
% LambdaX, LambdaY, kV, Xmove, Ymove
err = zeros(number,1);

for i = 1:number

% no wall effect
lambdaX = 1.209827074460249e+03; lambdaY = 1.857790332589524e+03;
kV = 3.9410; Xmove = 0.0253; Ymove = 0.1487; alpha = 0.6147;

% % with wall effect
% lambdaX = 1.173199043348060e+03; lambdaY = 1.709992650398775e+03;
% kV = 3.1353; Xmove = 0.4747; Ymove = 0.2259; alpha = 0.6145;

Sigma = [lambdaX,0; 0,lambdaY];
depx = Xc(i); depy = Yc(i); c = [depx, depy]; 
% the = atan((depy-Pe(2))/(depx-Pe(1)));
depH = H_data(:,:,i);    
V =  7.5889e+04 * Vol_data(i);

% ModelH = function_input_2d(X,Y,c-[Xmove,Ymove],kV*V,Sigma,the,xf,yr,yl);
ModelH = unit_sym_input_2d(X,Y,c-[Xmove,Ymove],kV*V,Sigma,alpha);

% figure
% subplot(1,2,1)
% mesh(X,Y,ModelH)
% xlim([0 170])
% zlim([-5 40])
% 
% % figure
% subplot(1,2,2)
% mesh(X,Y,depH)
% xlim([0 170])
% zlim([-5 40])

%% Error analysis
error_before = immse(depH, ModelH);
err(i) = sqrt(error_before);
% err(i)

%% Optimization
theta0 = [2000,2000,3,0,0,1];
lb = [1000,1000,2.5,0,-10,0.6];
ub = [3000,3000,4,20,10,2.0];

% fun = @(theta)immse(depH, function_input_2d(X,Y,c-[theta(4),theta(5)],theta(3)*V,[theta(1),0;0,theta(2)],the,xf,yr,yl));
fun = @(theta)immse(depH, unit_sym_input_2d(X,Y,c-[theta(4),theta(5)],theta(3)*V,[theta(1),0;0,theta(2)],theta(6)));

options = optimset('Display','iter','PlotFcns',@optimplotfval);

tic
theta = fmincon(fun,theta0,[],[],[],[],lb,ub);
% theta = fminsearch(fun,theta0);
toc

parameters(:,i) = theta';

%% result
lambdaX = theta(1); lambdaY = theta(2); kV = theta(3);
Xmove = theta(4); Ymove = theta(5); alpha = theta(6);

Sigma = [lambdaX,0; 0,lambdaY];
% ModelH = function_input_2d(X,Y,c-[Xmove,Ymove],kV*V,Sigma,the,xf,yr,yl);
ModelH = unit_sym_input_2d(X,Y,c-[Xmove,Ymove],kV*V,Sigma,alpha);

% error_before
error_after = immse(depH, ModelH);
error_before - error_after;

% figure
% % subplot(1,2,1)
% mesh(X,Y,ModelH + H0)
% xlim([0 170])
% zlim([-50 10])
% % zlim([-5 35])
% 
% figure
% % subplot(1,2,2)
% mesh(X,Y,depH + H0)
% xlim([0 170])
% zlim([-50 10])
% % zlim([-5 35])

% 
% input('Next data?[Enter]')
% close all
end

%% summary
maybeGoodResult = mean(parameters,2)
% resultant parameters are obtained from 7 data

LLambdaX = parameters(1,:);
LLambdaY = parameters(2,:);
kkV = parameters(3,:);
XXmove = parameters(4,:);
YYmove = parameters(5,:);
AAlpha = parameters(6,:);

lambdaX = maybeGoodResult(1); lambdaY = maybeGoodResult(2);
kV = maybeGoodResult(3); 
Xmove = maybeGoodResult(4); Ymove = maybeGoodResult(5);
alpha = maybeGoodResult(6);
% save nominal_modle_para_centerMove lambdaX lambdaY kV Xmove Ymove

mean(err)
figure; hold on; grid on;
plot(err,'LineStyle','none','Marker','x','Color','[0, 0, 0]','LineWidth',1.2);
plot([0 number],[mean(err) mean(err)],'--','Linewidth',1.2)
xlim([0 number]);   ylim([0 8]);
xlabel("Data No.")
ylabel("nMSE[mm]");