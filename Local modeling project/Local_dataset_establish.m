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
Xc = [Xc;62.5;37.5;25;25;37.5;62.5;75;75;50;50;62.5;37.5;37.5;37.5;62.5;62.5]*1.7;
Yc = [Yc;40;40;20;-20;-40;-40;-20;20;20;-20;0;0;-20;20;20;-20];

%% 3. The error model

% the center, volume -> error data
dep_center = [Xc,Yc]; Vol = zeros(number,1);
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

%% 4.1 Determine the local area for each data
% Output:
% H_local(m,n,number) - the data from H_error, outside local area is all 0.
% X_local(m,n,number) - the data of X, outside local area is all 0.
% Y_local(m,n,number) - the data of Y, outside local area is all 0.

% Parameters
ux = 70; uy = 60;   % area upper bound to center
lx = 50; ly = 60;   % area lower bound to center

% Get the area of local space
X_range = [Xc-lx, Xc+ux]; Y_range = [Yc-ly, Yc+uy];
Local_area = ones(m,n,number); % 1 and 0 meshgrid, 1 represents the local area
H_local = zeros(m,n,number); % The data of local area
X_local = zeros(m,n,number);
Y_local = zeros(m,n,number); % the meshgrid of position value for local area

for k = 1:number
    % represent with logical value
    X_area = X(1,:)<X_range(k,1); X_area = X_area | X(1,:)>X_range(k,2);
    Y_area = Y(:,1)<Y_range(k,1); Y_area = Y_area | Y(:,1)>Y_range(k,2);
    
    % convert logical value to double
    for i = 1:length(Y_area)
        for j = 1:length(X_area)
            if X_area(j) == 0 && Y_area(i) == 0
                Local_area(i,j,k) = 1;
            else
                Local_area(i,j,k) = 0;
            end
        end
    end
    
    % extract for X and Y
    X_local(:,:,k) = Local_area(:,:,k) .* X; 
    Y_local(:,:,k) = Local_area(:,:,k) .* Y;
    % extract the local data from data, not local => 0
    H_local(:,:,k) = Local_area(:,:,k) .* H_error(:,:,k);
end

%% 5. Make it sparse!
% Because different data has different size of local area,
% here let us make it sparse first.

step_length = 15; % level of sparseness

H_local_sparse = []; X_local_sparse = []; Y_local_sparse = [];

ix = 1; jy = 1;
for i = 1:step_length:m
    for j = 1:step_length:n
        H_local_sparse(ix,jy,:) = H_local(i,j,:);
        X_local_sparse(ix,jy,:) = X_local(i,j,:);
        Y_local_sparse(ix,jy,:) = Y_local(i,j,:);
        jy = jy + 1;
    end
    ix = ix + 1;
    jy = 1;
end

%% 6. Plot to check the local area and data
num = 2; % the data number: 1~61

% Evaluation: plot H_local to see whether your local area covers well!!
figure; mesh(H_error(:,:,num)); figure; mesh(H_local(:,:,num));

% Evaluation idea: Show the local area
figure; hold on;
title("Data No." + num + ", the dataset point(red) and selected point(black)")

% full data points
plot(X,Y,'LineStyle','none','Marker','.','Color','[0.8500, 0.3250, 0.0980]')

% soil loading center
plot(Xc(num),Yc(num),'LineStyle','none','Marker','x','Color','[0, 0, 0]','LineWidth',1.7)

% local area bounds
plot([X_range(num,1),X_range(num,2)],[Y_range(num,1),Y_range(num,1)],'Color','[0, 0, 0]')
plot([X_range(num,1),X_range(num,2)],[Y_range(num,2),Y_range(num,2)],'Color','[0, 0, 0]')
plot([X_range(num,1),X_range(num,1)],[Y_range(num,1),Y_range(num,2)],'Color','[0, 0, 0]')
plot([X_range(num,2),X_range(num,2)],[Y_range(num,1),Y_range(num,2)],'Color','[0, 0, 0]')

% selected data points
plot(X_local_sparse(:,:,num), Y_local_sparse(:,:,num),'LineStyle','none','Marker','x','Color','[0, 0, 0]')
hold off; % (0,0) are not aovided.

%% 7. Make a dataset
% pick up all the selected data points in the determined local area!
% and also their normalized relative positions.
% X_data = [Xc, Yc, Xre, Yre]
% Y_data = H

[m,n,number] = size(H_local_sparse);
X_data = [];
Y_data = [];

for k = 1:number
    
    % get rid of the 0 elements in H_local_sparse, X_local_sparse, Y_local_sparse
    XX = X_local_sparse(:,:,k);
    YY = Y_local_sparse(:,:,k);
    HH = H_local_sparse(:,:,k);
    for i = m:-1:1
        if sum(XX(i,:))==0 && sum(YY(i,:))==0 && sum(HH(i,:))==0
            XX(i,:) = []; YY(i,:) = []; HH(i,:) = [];
        end
    end
    for j = n:-1:1
        if sum(XX(:,j))==0 && sum(YY(:,j))==0 && sum(HH(:,j))==0
            XX(:,j) = []; YY(:,j) = []; HH(:,j) = [];
        end
    end
    
    % the relative positions and normalize
    Xre = (XX - Xc(k))/170; Yre = (YY - Yc(k))/80;
    Xre = Xre(:); Yre = Yre(:);
    Xc_nor = Xc(k)/170 * ones(length(Xre),1);
    Yc_nor = Yc(k)/80 * ones(length(Yre),1);
    
    % save the data to dataset vector
    Y_data = [Y_data; HH(:)];
    X_data = [X_data; Xc_nor, Yc_nor, Xre, Yre];
end

disp("The input data is N = " + length(X_data(:,1)) + " with D = " + length(X_data(1,:)))

% X_test
X_test = [0.5*ones(9400,1), zeros(9400,1), (X(:)-85)/170, Y(:)/80];

save local_dataset X_data Y_data X Y