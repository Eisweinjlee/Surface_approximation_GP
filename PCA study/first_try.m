%% Principle component analysis study
% We want to use PCA method to present the property of images (soil dis.)
% Author: Yang Li
% Date: Dec 20th, 2019
clear
close all

%% 1. Load some data
% excavator_data;

docName1 = "Approx_Surface\ErrorData training project\dataset_20191108_normalized\";
docName2 = "Approx_Surface\ErrorData training project\dataset_20191209_normalized\";

m = 94; n = 100;
H_data = zeros(m,n,97);

% load the empty vessel
load(docName1 + "00.mat",'dep')
H0 = dep;

% load docName1 - Position difference
number = 0;
% H_data(:,:,number) = H0;

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

%% 2. Vectorization
data_vector = zeros(9400,number);

for i = 1:number
    data_vector(:,i) = reshape(H_data(:,:,i),[9400,1]);
end

% mean value of each data
vmean = mean(data_vector,1);
M = data_vector - vmean; % vbar1, vbar2, ....

%% 3. Singular value decomposition (SVD)

[U,S,V] = svd(M);

%% 4. Approximation
err = zeros(number,1);
i = 1;
for r = 1:number
% r = 5; % see the S, we choose r = 5

cbar = S*V';
cbarre = cbar(1:r,:); % Principal component of each profile

Ure = U(:,1:r); % reduced weight matrix

M_approx = Ure * cbarre;
err(i) = immse(M,M_approx);
i=i+1;
end

plot(err)
xlabel("$$r'$$",'Interpreter','latex')
ylabel("Error")
grid on

%% 5. PCA
r = 4;  % is enough for our project to use r=4
cbar = S*V';

cbarre = cbar(1:r,:); % Principal component of each profile
Ure = U(:,1:r); % reduced weight matrix

%% 6. show the result

M_approx = Ure * cbarre;
H_approx = M_approx + vmean;

for i = [3 8 18]
    H_plot = reshape(H_approx(:,i),[94,100]);
    figure; subplot(1,2,1);mesh(H_data(:,:,i));
    subplot(1,2,2);mesh(H_plot);
end