%% Try the density of data
% Author: LI, Yang
% Date: Dec 2nd, 2019
close all
clear

%% 1. Load the 7 data
load('data_seven_centerdep/Error_data1011.mat')
H = H_error;
clear H_error
m = length(H(:,1,1));
n = length(H(1,:,1));
number = length(H(1,1,:));

H_full = H; % Full dataset


%% 2. Make the data local
load C:\Users\LI\Desktop\Approx_Surface\Vessel_XY.mat

Xc = 85; Yc = 0;

ux = 80; uy = 60;   % upper bound
lx = 40; ly = 60;   % lower bound

% Get the area of local space
X_range = [Xc - lx; Xc + ux]; Y_range = [Yc - ly; Yc + uy];
Local_area = ones(m,n);
% represent with logical value
X_area = X(1,:)<X_range(1); X_area = X_area | X(1,:)>X_range(2);
Y_area = Y(:,1)<Y_range(1); Y_area = Y_area | Y(:,1)>Y_range(2);
% convert logical value to double
for i = 1:length(Y_area)
    for j = 1:length(X_area)  
        if X_area(j) == 0 && Y_area(i) == 0
            Local_area(i,j) = 1;
        else
            Local_area(i,j) = 0;
        end
    end
end

% extract the local data from full data
for k = 1:number
AA = Local_area .* H_full(:,:,k);
BB = AA;

for i = length(Y_area):-1:1
   if  Y_area(i) == 1
       BB(i,:) = [];
   end
end
for i = length(X_area):-1:1
   if  X_area(i) == 1
       BB(:,i) = [];
   end
end

% figure
% mesh(BB)
% title("No."+k)

H_local(:,:,k) = BB;
end

% extract for X and Y
X_local = Local_area .* X;
Y_local = Local_area .* Y;

for i = length(Y_area):-1:1
   if  Y_area(i) == 1
       X_local(i,:) = [];
       Y_local(i,:) = [];
   end
end
for i = length(X_area):-1:1
   if  X_area(i) == 1
       X_local(:,i) = [];
       Y_local(:,i) = [];
   end
end

figure
title("The dataset point(red) and selected point(black)")
hold on
plot(X,Y,'LineStyle','none','Marker','.','Color','[0.8500, 0.3250, 0.0980]')
plot(Xc,Yc,'LineStyle','none','Marker','x','Color','[0, 0, 0]','LineWidth',1.7)

plot([X_range(1),X_range(2)],[Y_range(1),Y_range(1)],'Color','[0, 0, 0]')
plot([X_range(1),X_range(2)],[Y_range(2),Y_range(2)],'Color','[0, 0, 0]')
plot([X_range(1),X_range(1)],[Y_range(1),Y_range(2)],'Color','[0, 0, 0]')
plot([X_range(2),X_range(2)],[Y_range(1),Y_range(2)],'Color','[0, 0, 0]')

%% 3. Make it sparse
step_length = 8; % level of sparseness

H_local_sparse = [];
X_local_sparse = [];
Y_local_sparse = [];

[mm,nn,kk] = size(H_local);
ix = 1; jy = 1;

for i = 1:step_length:mm
    for j = 1:step_length:nn
        H_local_sparse(ix,jy,:) = H_local(i,j,:);
        X_local_sparse(ix,jy) = X_local(i,j);
        Y_local_sparse(ix,jy) = Y_local(i,j);
        jy = jy + 1;
    end
    ix = ix + 1;
    jy = 1;
end

plot(X_local_sparse,Y_local_sparse,'LineStyle','none','Marker','x','Color','[0, 0, 0]')

% figure
% mesh(X_local,Y_local,H_local(:,:,1))

%% 4. Make it a data set

% normalize the X = [0,1] and Y = [-1,1]
X_local_sparse = X_local_sparse./170;
Y_local_sparse = Y_local_sparse./80;

% Matrices to verctors
X_data_base = [X_local_sparse(:),Y_local_sparse(:)];
X_data = [];
for i = 1:number
   X_data = [X_data;X_data_base]; 
end

Y_data = H_local_sparse(:);

X_local = X_local./170;
Y_local = Y_local./80;
X_test = [X_local(:),Y_local(:)];

save local_sparse_7_data X_data Y_data X_test X_local Y_local
disp("The dataset size is N=" + length(X_data(:,1)) + " with D=" + length(X_data(:,2)))