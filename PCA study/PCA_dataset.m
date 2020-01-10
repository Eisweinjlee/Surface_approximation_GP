%% Dataset for PCA
% the dataset for the principal components calculation
% Author: Yang Li
% Date: Jan 10th, 2020
clear

%% 1. Loading the data
docName1 = "Approx_Surface\ErrorData training project\dataset_20191108_normalized\";
docName2 = "Approx_Surface\ErrorData training project\dataset_20191209_normalized\";
docName3 = "Approx_Surface\ErrorData training project\dataset_20200108_normalized\";

load(docName1 + "00.mat",'dep')
H0 = dep; [m,n] = size(dep);

number = 0;
for i = 1:9
    number = number + 1;
    filename = docName1 + num2str(10*i+3)+".mat";
    load(filename)
    H_data(:,:,number) = dep;
end
for i = 1:16
    number = number + 1;
    filename = docName1 + num2str(100+i)+".mat";
    load(filename)
    H_data(:,:,number) = dep;
end
for i = 1:9
    for j = [1 2 4]
        number = number + 1;
        filename = docName2 + num2str(10*i+j)+".mat";
        load(filename)
        H_data(:,:,number) = dep;
    end
end
% b10 ~ b61
for i = 1:6
    filename = docName3 + "b" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% bb10 ~ bb81
for i = 1:8
    filename = docName3 + "bb" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% c10 ~ c91
for i = 1:9
    filename = docName3 + "c" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% l10 ~ l61
for i = 1:6
    filename = docName3 + "l" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% r10 ~ r61
for i = 1:6
    filename = docName3 + "r" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% t10 ~ t61
for i = 1:6
    filename = docName3 + "t" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% tt10 ~ tt81
for i = 1:8
    filename = docName3 + "tt" + num2str(10*i+1)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    number = number + 1;
end
% 11 ~ 103
for i = 1:10
    filename = docName3 + num2str(10*i+2)+".mat"; load(filename);
    H_data(:,:,number + 1) = dep;
    filename = docName3 + num2str(10*i+3)+".mat"; load(filename);
    H_data(:,:,number + 2) = dep;
    number = number + 2;
end

%% 2. Extract the small area
x1 = (2:24); x2 = (25:47); x3 = (48:70); x4 = (71:93);
y1 = (1:25); y2 = (26:50); y3 = (51:75); y4 = (76:100);

H_small = zeros(46,50,number*9);
num_small = 0;
for i = 1:number
    H_small(:,:,num_small+1) = H_data([x1,x2],[y1,y2],i);
    H_small(:,:,num_small+2) = H_data([x1,x2],[y2,y3],i);
    H_small(:,:,num_small+3) = H_data([x1,x2],[y3,y4],i);
    H_small(:,:,num_small+4) = H_data([x2,x3],[y1,y2],i);
    H_small(:,:,num_small+5) = H_data([x2,x3],[y2,y3],i);
    H_small(:,:,num_small+6) = H_data([x2,x3],[y3,y4],i);
    H_small(:,:,num_small+7) = H_data([x3,x4],[y1,y2],i);
    H_small(:,:,num_small+8) = H_data([x3,x4],[y2,y3],i);
    H_small(:,:,num_small+9) = H_data([x3,x4],[y3,y4],i);
    num_small = num_small + 9;
end

%% 3. Vectorization
data_vector = zeros(46*50,num_small);

for i = 1:num_small
    data_vector(:,i) = reshape(H_small(:,:,i),[46*50,1]);
end

% mean value of each data
vmean = mean(data_vector,1);
M = data_vector - vmean; % vbar1, vbar2, ....

%% 3. Singular value decomposition (SVD)
tic
[U,S,V] = svd(M);
toc

%% 4. Find a good r'
err = zeros(10,1);
i = 1;
for r = 1:10

cbar = S*V';
cbarre = cbar(1:r,:); % Principal component of each profile

Ure = U(:,1:r); % reduced weight matrix

M_approx = Ure * cbarre;
err(i) = immse(M,M_approx);
i = i+1;
end

plot(err)
xlabel("$$r'$$",'Interpreter','latex')
ylabel("RMSE")
grid on

%% 5. PCA
r = 4;  % is enough for our project to use r = 3
cbar = S*V';

cbarre = cbar(1:r,:); % Principal component of each profile
Ure = U(:,1:r); % reduced weight matrix

M_approx = Ure * cbarre;
H_approx = M_approx + vmean;

% for i = [198 584 720 900 1023]
%     H_plot = reshape(H_approx(:,i),[46,50]);
%     figure; subplot(1,2,1);mesh(H_small(:,:,i));
%     subplot(1,2,2);mesh(H_plot);
% end

% for i = 1:num_small
%     H_plot = reshape(H_approx(:,i),[46,50]);
%     rmse(i) = sqrt(immse(H_small(:,:,i),H_plot));
% end
% figure;plot(rmse)