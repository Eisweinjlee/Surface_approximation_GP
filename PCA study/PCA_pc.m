function [H_small_PC, H_small_mean] = PCA_pc(Xc,H)
%% Principal components calculation for the area around loading center.
% The dataset is obtained to be M with PCA_dataset.m
% Author: Yang Li
% Date: Jan 11th, 2020

% Input: 
% Xc - the loading center - 2x1 vector
% H  - the shape before loading - 94x100

% % for testing
% clear
% load("H_test")
% Xc = [80;-20];

%% 1. Load the PCA_dataset M
load("Approx_Surface\PCA study\M_matrix_PCA.mat")

%% 2. pick the small area around loading center
depx = round(Xc(1)*100/170); depy = round((Xc(2)+80)*94/160);
x_range = ((depy-22):(depy+23)); y_range = ((depx-24):(depx+25));
H_small = H(x_range,y_range);

%% 3. Vectrerize and put into M & mean value
H_small_vec = H_small(:);
H_small_mean = mean(H_small_vec);

v_bar = H_small_vec - H_small_mean;

M = [M, v_bar];

%% 4. Singular value decomposition (SVD)
[U,S,V] = svd(M);

%% 5. PC
r = 4;  % is enough for our project to use r = 4
cbar = S*V';

% cbarre = cbar(1:r,:); % Principal component of each profile
% Ure = U(:,1:r); % reduced weight matrix

H_small_PC = cbar(1:r,end);
end