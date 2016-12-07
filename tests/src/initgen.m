% Generate initializer matrices for the hierarchical clustering
% code.

num_matrices = 64;

% Reuters matrix: 12411 x 7984
% W matrices are 12411 x 2
% H matrices are 2 x 7984
%m = 12411;
%n =  7984;

% 20 news matrix: 39727 x 11237
% W matrices are 39727 x 2
% H matrices are 2 x 11237
m = 39727
n = 11237

% Wikipedia big input matrix: 162407 x 85762
% W matrices are 162407 x 2
% H matrices are 2 x 85762
%m = 162407
%n = 85762

%if ~exist('matrices.reuters', 'dir')
%    mkdir('matrices.reuters');
%end

if ~exist('matrices.20news', 'dir')
    mkdir('matrices.20news');
end

%if ~exist('matrices.wikipedia', 'dir')
%    mkdir('matrices.wikipedia');
%end


for i = 1:num_matrices
    W = rand(m, 2);
    H = rand(2, n);
    
    %w_file = sprintf('matrices.wikipedia/Winit_%d.csv', i);
    w_file = sprintf('matrices.20news/Winit_%d.csv', i);
    %w_file = sprintf('matrices.reuters/Winit_%d.csv', i);
    csvwrite(w_file, W);
    
    %h_file = sprintf('matrices.wikipedia/Hinit_%d.csv', i);
    h_file = sprintf('matrices.20news/Hinit_%d.csv', i);
    %h_file = sprintf('matrices.reuters/Hinit_%d.csv', i);
    csvwrite(h_file, H);
end
