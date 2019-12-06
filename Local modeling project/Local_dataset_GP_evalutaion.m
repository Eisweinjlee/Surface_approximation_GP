%% Local_dataset_GP_evaluation
% The trained GP model for error distribution is evaluated.
% Author: Li, Yang
% Date: Dec 5th, 2019
close all
clear

%% 1. Load the model and specify parameters

% load the dataset: X_data, Y_data, X, Y
load local_dataset.mat
% load trained parameters: hyp_sparseGP
load data20191204.mat

% Specify the mean, cov, likelihood
meanfunc = [];                      % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;           % inference with Guassian Likelihood

% inducing points
xu = X_data(1:10:end,:); cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s', 0));
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1

%% 2. The prediction

% normalized center: Xc=[0,1], Yc=[-1,1]
Xc_nor = 0.25; Yc_nor = -0.5;
Xc = Xc_nor * 170; Yc = Yc_nor * 80;
% test data
X_test = [Xc_nor*ones(9400,1), Yc_nor*zeros(9400,1), (X(:)-Xc)/170, (Y(:)-Yc)/80];

% prediction
[ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
    X_data, Y_data, X_test);

% reshape the result
m = 94; n = 100;
H_error_pred = reshape(ymu,[m,n]);
CovDist = reshape(ys2,[m,n]);

% % Plot: the mean prediction
% figure; mesh(X, Y, H_error_pred)
% title("$X_c=$"+Xc*170+", $Y_c=$"+Yc*80,'Interpreter','latex')

% % Plot: covariance analysis
% s_upper = ymu + 2*sqrt(ys2); s_lower = ymu - 2*sqrt(ys2);
% s_upper = reshape(s_upper,[m,n]); s_lower= reshape(s_lower,[m,n]);
% figure; hold on;
% sur1 = surf(X,Y,H_error); sur1.EdgeColor = 'flat';
% sur2 = surf(X,Y,s_upper,'FaceAlpha',0.3); sur2.EdgeColor = 'none';
% sur3 = surf(X,Y,s_lower,'FaceAlpha',0.3); sur3.EdgeColor = 'none';
% view([-23,20]); hold off;
%
% figure; mesh(X,Y,CovDist);

%% 3. Relative covariance modification

% find the data close to center X(1,:), Y(:,1)
disX = (X(1,:)-Xc).^2; disY = (Y(:,1)-Yc)'.^2;
[minX_v,minX_p] = min(disX); [minY_v,minY_p] = min(disY);
data_cen = [X(1,minX_p);Y(minY_p,1)];

% % Plot: the loading center and data center
% figure;hold on;
% plot(X,Y,'LineStyle','none','Marker','.','Color',[218,165,32]/255)
% plot(Xc*170,Yc*80,'LineStyle','none','Marker','+','Color',[199,21,133]/255,'LineWidth',1.5)
% plot(data_cen(1),data_cen(2),'LineStyle','none','Marker','+','Color',[0,0,128]/255,'LineWidth',1.5)
% hold off;

% calculate all relative covariance
cov_cen = CovDist(minY_p,minX_p); % the covariance of center
w = cov_cen ./ CovDist;

% Sigmoid function
a = 2; b = 7; % https://www.desmos.com/calculator/kn9tpwdan5
% a = 9; b = 14.2;
gw = 1./(1 + exp(a + -b.*w)); % sigmoid
% q = 0:0.01:max(w,[],'all'); % plot to see the sigmoid
% gw = 1./(1 + exp(a + -b.*(q-c))); % sigmoid
% figure; plot(q,gw); grid on

% modify the mean predition through sigmoid
H_modified = gw .* H_error_pred;
% % Plot: before modification and after
% figure; subplot(1,2,1); mesh(X,Y,H_error_pred); zlim([-20 20]); title("before");
% subplot(1,2,2); mesh(X,Y,H_modified); zlim([-20 20]);title("after");

%% 4. Add to the nominal Gaussian pdf model
excavator_data;
Vol = 7.5889e+04;

H_nominal = function_input_2d(X,Y,[Xc;Yc],3.75*Vol,Sigma,the,xf,yr,yl);
H_noGP = H0 + H_nominal;
H_combination = H_noGP + H_modified;

%% 5. Evaluate the performance
docName = "Approx_Surface\ErrorData training project\dataset_20191108_normalized\";
filename = docName + "93" +".mat";
load(filename)
err = sqrt(immse(H_combination, dep))

figure; 
subplot(1,3,1); mesh(X,Y,H_noGP); 
xlim([0 170]); zlim([-50 40]); title("Nominal");zlabel("h[mm]")
subplot(1,3,2); mesh(X,Y,H_combination); 
xlim([0 170]); zlim([-50 40]); title("Nominal+GP");
subplot(1,3,3); mesh(X,Y,dep); zlim([-20 20]);
xlim([0 170]); zlim([-50 40]); title("Real, nMSE= "+err);

% H_data = zeros(m,n,61);

% % load others
% number = 0;
% for i = 1:9
%     for j = 1:5
%         number = number + 1;
%         filename = docName + num2str(10*i+j)+".mat";
%         load(filename)
%         H_data(:,:,number) = dep - H0;
%     end
% end
% for i = 1:16
%     number = number + 1;
%     filename = docName + num2str(100+i)+".mat";
%     load(filename)
%     H_data(:,:,number) = dep - H0;
