% *Remove specified figure:*
%%
% _Deletes the current figure or the specified figure(s)._

close all

%%
% *Remove items from workspace, freeing up system memory:*
%%
% _Removes all variables from the current workspace, releasing them from_
% _system memory._

clear all

%%
% *Clear Command Window:*
%%
% _Clears all input and output from the Command Window display, giving you_
% _a "cleanscreen"._

clc
% Set the number of nodes and region size
n=input(['Enter the number of nodes you assumed,\n', ...
    'the value should be integer more than 1 \n']);  
region_size = 30;

% Generate random coordinates for the nodes
x = rand(1, n) * region_size;
y = rand(1, n) * region_size;

% Generate random traffic demands for the requests
m_lambda = 150; % Average traffic amount per request
requests = cell(n, 1);

for i = 1:n
    k = poissrnd(1); % Number of requests originating from node i
    destinations = setdiff(1:n, i); % Pick destinations from other nodes
    dest_indices = randsample(destinations, k, true);
    traffic_demands = normrnd(m_lambda, sqrt(0.5 * m_lambda), k, 1);
    
    requests{i} = [i * ones(k, 1), dest_indices', traffic_demands];
end

% Plot the nodes and requests
scatter(x, y, 'filled');
hold on;
for i = 1:n
    for j = 1:size(requests{i}, 1)
        s = requests{i}(j, 1);
        d = requests{i}(j, 2);
        plot([x(s), x(d)], [y(s), y(d)], 'b');
    end
end
hold off;
axis([0 region_size 0 region_size]);
title('Topology for Non-Splittable and Splittable Traffic');
xlabel('X');
ylabel('Y');