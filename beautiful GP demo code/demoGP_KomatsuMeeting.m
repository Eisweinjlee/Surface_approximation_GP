%% Demo of Gaussian process for Komatsu meeting
% LI Yang// Nov 13th, 2019
clear
close all

% Method selection: (0)VFE, (1)FTIC
selection_flag = 0;

% New data(0) or load data(1)
data_type_flag = 0;


%% The data generation

% latent function
x = linspace(-2,2,100);
y = 0.6*sin(2*pi.*x + 0.2) + 0.8*cos (pi.*x -0.5);
% plot(x,y)
% axis([-2 2 -2 2])

[m,n] = size(y);

% Gaussian noise
var = 0.5;
y1 = y + var .* randn(m,n);
y2 = y + var .* randn(m,n);
y3 = y + var .* randn(m,n);
y4 = y + var .* randn(m,n);

% Cauchy noise
cau = 0.02;
y5 = y + cau * tan(pi*(rand(m,n)-1/2));

%% Dataset

X_data_full = [x';x';x';x';x'];
Y_data_full = [y1';y2';y3';y4';y5'];

numberOfIndices = 215;
random_indices = randi(500,[1,numberOfIndices]);

% plot(X_data,Y_data,'*')
% axis([-2 2 -2 2])

X_test = x';

for i = 5:70:numberOfIndices
% for i = 215
    
    X_data = X_data_full(random_indices(1:i));
    Y_data = Y_data_full(random_indices(1:i));
    
    %% Gaussian Process
    
    N = 200;  % iteration times limitation
    
    % Specify the mean, cov, likelihood
%     meanfunc = [];                    % empty: don't use a mean function
    meanfunc = {@meanSum, {@meanLinear, @meanConst}}; 
    hyp_init.mean = [0; 0.5];
    
    covfunc = {@covSEard};              % ARD SE
%     covfunc = {@covMaternard,3};        % Matern ARD (non-smooth)
    hyp_init.cov = [0; 0];

    likfunc = {@likGauss};              % Gaussian likelihood
    infmethod = @infGaussLik;               % inference with Guassian Likelihood
    hyp_init.lik = -1;
    
    % The hyperparameter struct
%     hyp_init = struct('mean', [], 'cov', [0 0], 'lik', -1);
    
    % Optimization hyperparameters
    tic;
    hyp_exactGP = minimize(hyp_init, @gp, -N, infmethod, meanfunc, covfunc,...
        likfunc, X_data, Y_data);
    Time_of_ExactGP = toc;
    
    % The dense prediction
    [mu, s2] = gp(hyp_exactGP, infmethod, meanfunc, covfunc, likfunc,...
        X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest
    
    %
    % % inducing points
    % % mid = floor(length(X_test)/2);
    % % xu = X_test(mid-3:mid+3); % (1) bad choice of initial inducing set
    % xu = X_test(1:30:end); % (2) generous choice of initial inducing set
    % cov = {'apxSparse', covfunc, xu};
    % inff = @(varargin) infmethod(varargin{:},struct('s', selection_flag));
    % % VFE, opt.s = 0; SPEP, 0 <opt.s < 1; FITC, opt.s = 1
    %
    % tic; hyp_init.xu = xu;
    % hyp_sparseGP = minimize(hyp_init, @gp, -N, inff, meanfunc, cov, likfunc,...
    %     X_data, Y_data);
    % Time_of_SparseGP = toc;
    %
    % [ymu,ys2] = gp(hyp_sparseGP, inff, meanfunc, cov, likfunc,...
    %     X_data, Y_data, X_test);
    %
    % % % Marginal likelihood and derivatives
    % % [nlZ,dnlZ] = gp(hyp_sparseGP,infmethod, meanfunc, cov,...
    % %     likfunc, X_data, Y_data);
    
    
    
    %% Plot
    
    % Analysis
    Time_of_ExactGP
    % Time_of_SparseGP
    
    exactGP_mse = immse(mu,y')
    f1 = [mu+2*sqrt(s2), flipdim(mu-2*sqrt(s2),1)];
    
    % sparseGP_mse = immse(ymu,y')
    % f2 = [ymu+2*sqrt(ys2), flipdim(ymu-2*sqrt(ys2),1)];
    
    figure
    hold on
    % Dataset
    plot(X_data,Y_data,'rx','Linewidth',1)
    plot(x,y,'g','Linewidth',1.5)
    
    % Predictive mean plot
    plot(X_test,mu,'b','Linewidth',1.5)
    % plot(X_test,ymu,'r','Linewidth',1.5)
    
    % exact GP predictive variance
    % fill([X_test; flipdim(X_test,1)], f1(:), [0.9 0.9 0.9]);
    plot(X_test,f1(:,1),'--','Linewidth',1,'Color','#0072BD')
    plot(flipdim(X_test,1),f1(:,2),'--','Linewidth',1,'Color','#0072BD')
    
    titlename = "Number of datapoint: " + i;
    title(titlename)
    
%     % sparse GP predictive variance
%     % fill([X_test; flipdim(X_test,1)], f2(:), [0.9 0.9 0.9]);
%     plot(X_test,f2(:,1),'--','Linewidth',1,'Color','#A2142F')
%     plot(flipdim(X_test,1),f2(:,2),'--','Linewidth',1,'Color','#A2142F')
%     
%     % Inducing point displacement
%     plot(xu,2.5*ones(length(xu),1),'r^') % initial
%     plot(hyp_sparseGP.xu,-2.5*ones(length(hyp_sparseGP.xu),1),'rv') % result
    
    hold off
    grid on
    axis([-2 2 -2.5 2.5])
    
    % figure
    % K = feval(covfunc{:}, hyp_exactGP.cov, sort(X_data));
    % contour(K,'Fill','on')
end

%% Plot out_of_dataset
X_test = linspace(-3.5,3.5,500)';
% The dense prediction
[mu, s2] = gp(hyp_exactGP, infmethod, meanfunc, covfunc, likfunc,...
    X_data, Y_data, X_test);  % Xtrain,Ytrain,Xtest

f1 = [mu+2*sqrt(s2), flipdim(mu-2*sqrt(s2),1)];

figure
hold on
% Dataset
plot(X_data,Y_data,'rx','Linewidth',1)
plot(x,y,'g','Linewidth',1.5)

% Predictive mean plot
plot(X_test,mu,'b','Linewidth',1.5)
% plot(X_test,ymu,'r','Linewidth',1.5)

% exact GP predictive variance
% fill([X_test; flipdim(X_test,1)], f1(:), [0.9 0.9 0.9]);
plot(X_test,f1(:,1),'--','Linewidth',1,'Color','#0072BD')
plot(flipdim(X_test,1),f1(:,2),'--','Linewidth',1,'Color','#0072BD')

% % sparse GP predictive variance
% % fill([X_test; flipdim(X_test,1)], f2(:), [0.9 0.9 0.9]);
% plot(X_test,f2(:,1),'--','Linewidth',1,'Color','#A2142F')
% plot(flipdim(X_test,1),f2(:,2),'--','Linewidth',1,'Color','#A2142F')

% % Inducing point displacement
% plot(xu,2.5*ones(length(xu),1),'r^') % initial
% plot(hyp_sparseGP.xu,-2.5*ones(length(hyp_sparseGP.xu),1),'rv') % result

hold off
grid on
axis([-3.5 3.5 -2.5 2.5])
