% Establishing local dataset for deposit data
% with the "dataset_20191108_normalized"
% Author: Li, Yang
% Date: Dec 3rd, 2019

close all
clear

excavator_data;
docName = "Approx_Surface\ErrorData training project\dataset_20191108_normalized\";

%% 1. Loading the data (changed for different data)

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

%% 2. The loading center data (normalized)
Xc = [25*ones(5,1);50*ones(5,1);75*ones(15,1);50*ones(10,1);25*ones(10,1)];
Yc = [zeros(15,1);40*ones(5,1);-40*ones(10,1);40*ones(10,1);-40*ones(5,1)];
Xc = [Xc;62.5;37.5;25;25;37.5;62.5;75;75;50;50;62.5;37.5;37.5;37.5;62.5;62.5]./100;
Yc = [Yc;40;40;20;-20;-40;-40;-20;20;20;-20;0;0;-20;20;20;-20]./80;

%% 3. The error model

% the center, volume -> error data
dep_center = [Xc*170,Yc*80]; Vol = zeros(number,1);
H_error = zeros(m,n,number);
for i = 1:number
    Vol(i) = sum(H_data(:,:,i) ,'all');         % volume
    
    delta_H = function_input_2d(X,Y,dep_center(i,:)',2.96*Vol(i),Sigma,the,xf,yr,yl);
    H_error(:,:,i) = H_data(:,:,i) - delta_H;   % error data
    
%     close all % test use (study about the error distribution!)
%     figure
%     mesh(delta_H)
%     title("Gaussian pdf")
%     figure
%     mesh(H_error(:,:,i))
%     title("Error Model")
    
end

%% 4. Make it local to the center
% local model, also the coordinates should be local
% consider a way to change the local area relative to the center


%% 5. Make it sparse!
