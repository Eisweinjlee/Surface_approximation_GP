%% Local_dataset_GP_ver1
% Author: Li, Yang
% Date: Dec 4th, 2019

close all
clear

Mode_flag = 0; % 0: training & predict, 1: only predict

%% load the dataset
load local_dataset-12-Dec-2019.mat
% X_data, Y_data, X, Y

% test data
xx = 0.5; yy = 0; vol = 1.0; % normalized loading center
X_test = [xx*ones(9400,1), yy*zeros(9400,1), (X(:)-xx*170)/170,...
    (Y(:)-yy*80)/80, vol*ones(9400,1)];

%% Gaussian Process parameters

N = 400;  % iteration times limitation

% Specify the mean, cov, likelihood
meanfunc = [];                    % empty: don't use a mean function
covfunc = {@covSEard};              % ARD SE
likfunc = {@likGauss};              % Gaussian likelihood
infmethod = @infGaussLik;               % inference with Guassian Likelihood
% The hyperparameter struct
hyp_init = struct('mean', [], 'cov', zeros(1,length(X_data(1,:))+1), 'lik', -1);

% inducing points
require_num = 1.5*log10(length(X_data(:,1)))^(length(X_data(1,:)));
step_length = floor(length(X_data(:,1))/require_num);
xu = X_data(1:step_length:end,:);
cov = {'apxSparse', covfunc, xu};
inff = @(varargin) infmethod(varargin{:},struct('s', 0));
% VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1


if Mode_flag == 0 % training & predict
    %% Sparse Gaussian Process training
    tic; hyp_init.xu = xu;
    hyp_sparseGP = minimize(hyp_init, @gp, -N, inff, meanfunc, cov, likfunc,...
        X_data, Y_data);
    Time_of_SparseGP = toc
    
    save("model-"+date, 'hyp_sparseGP')
    disp("Trained model is saved.")
    
    %% Predict
    [ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
        X_data, Y_data, X_test);
    
    % % Marginal likelihood and derivatives
    % [nlZ,dnlZ] = gp(hyp_sparseGP,infmethod, meanfunc, cov,...
    %     likfunc, X_data, Y_data);
    
    % reshape the result
    m = 94; n = 100;
    H_local_sparse = reshape(ymu,[m,n]);
    
    % figure
    figure
    mesh(X, Y, H_local_sparse)
    
elseif Mode_flag == 1 % only predict
    %% Predict
    
    % load trained parameters
    load("model-11-Dec-2019.mat")
    
    for xx = [0.2 0.8]
        for yy = [-0.5 0.0 0.5]
            for vol = [0.3 0.6]
            
            % test data
            % xx = 0.1; yy = 0.8; % normalized loading center
            X_test = [xx*ones(9400,1), yy*zeros(9400,1), (X(:)-xx*170)/170,...
                (Y(:)-yy*80)/80, vol*ones(9400,1)];
            
            [ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
                X_data, Y_data, X_test);
            
            % reshape the result
            m = 94; n = 100;
            H_local_sparse = reshape(ymu,[m,n]);
            Cov_dist = reshape(ys2,[m,n]);
            
            % figure
            figure
            mesh(X, Y, H_local_sparse)
            title("$X_c=$"+xx*170+", $Y_c=$"+yy*80+", $V=$"+vol, 'Interpreter','latex')
            
            figure
            mesh(X, Y, Cov_dist)
            
%             s_upper = ymu + 2*sqrt(ys2);
%             s_lower = ymu - 2*sqrt(ys2);
%             s_upper = reshape(s_upper,[m,n]);
%             s_lower = reshape(s_lower,[m,n]);
%             
%             figure
%             hold on
%             sur1 = surf(X,Y,H_local_sparse);
%             sur1.EdgeColor = 'flat';
%             sur2 = surf(X,Y,s_upper,'FaceAlpha',0.3);
%             sur2.EdgeColor = 'none';
%             sur3 = surf(X,Y,s_lower,'FaceAlpha',0.3);
%             sur3.EdgeColor = 'none';
%             view([-23,20])
%             hold off;
            end
        end
    end
    
end