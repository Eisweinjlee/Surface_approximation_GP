function [X_data, Y_data] = datasetReduction(X, Y, H)

% The easiest idea to reduce the density of the data points
% X, Y, H are m-by-n

m = length(H(:,1));
n = length(H(1,:));
steps = 2;  % 1/4 of data

x = X(1:steps:m,1:steps:n);
y = Y(1:steps:m,1:steps:n);
h = H(1:steps:m,1:steps:n);

X_data = [x(:), y(:)];
Y_data = h(:);

end
