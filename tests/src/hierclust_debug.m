
% Add the folder xdata3/matlab to the matlabpath prior to running
% this file.  Do this with the command:
%
%    addpath path/to/xdata3/matlab

% number of clusters
k = 9

% A = [0 4 8  12 16 20;
%      1 5 9  13 17 21;
%      2 6 10 14 18 22;
%      3 7 11 15 19 23]

% W = [0 4  8 12;
%      1 5  9 13;
%      2 6 10 14;
%      3 7 11 15]

% H = [0 4 8  12 16 20;
%      1 5 9  13 17 21;
%      2 6 10 14 18 22;
%      3 7 11 15 19 23]

% dense test
% [A] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
%                'rnd_256_256.csv']);

% [W] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
%                'w_init_256_16.csv']);
% [H] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
%                'h_init_16_256.csv']);

% sparse test

% Reuters matrix
%[A] = mmread('/Users/richardboyd/repos/bitbucket/xdata_data/reuters.mtx');
%[Winit] = csvread('/Users/richardboyd/repos/bitbucket/xdata3/tests/src/matrices.reuters/Winit_1.csv');
%[Hinit] = csvread('/Users/richardboyd/repos/bitbucket/xdata3/tests/src/matrices.reuters/Hinit_1.csv');

% 20 news matrix
[A] = mmread('/Users/richardboyd/repos/bitbucket/gatech_share/NMF_20news_input_tf_reduced/reduced_matrix.mtx');
[Winit] = csvread('/Users/richardboyd/repos/bitbucket/xdata3/tests/src/matrices.20news/Winit_1.csv');
[Hinit] = csvread('/Users/richardboyd/repos/bitbucket/xdata3/tests/src/matrices.20news/Hinit_1.csv');

% Wikipedia big input matrix
%[A] = mmread('/Users/richardboyd/repos/bitbucket/gatech_share/NMF_wikipedia_big_input_reduced/reduced_matrix.mtx');
%[Winit] = csvread('/Users/richardboyd/repos/bitbucket/xdata3/tests/src/matrices.wikipedia/Winit_1.csv');
%[Hinit] = csvread('/Users/richardboyd/repos/bitbucket/xdata3/tests/src/matrices.wikipedia/Hinit_1.csv');

fprintf('A: %d x %d\n', size(A, 1), size(A, 2));
fprintf('Norm of Winit: %f\n', norm(Winit, 'fro'));
fprintf('Norm of Hinit: %f\n', norm(Hinit, 'fro'));

if (size(Hinit, 1) ~= size(Winit, 2))
    fprintf('Winit and Hinit do not agree on the value of the inner dim.\n');
    exit(-1);
end

trial_allowance = 3;
unbalanced = 0.1;
vec_norm = 2.0;
normW = true;
anls_alg = @anls_entry_rank2_precompute;
tol = 1e-4;
maxiter = 10000;

params = [];
params.trial_allowance = trial_allowance;
params.unbalanced = unbalanced;
params.vec_norm = vec_norm;
params.normW = normW;
params.anls_alg = anls_alg;
params.tol = tol;
params.maxiter = maxiter;

% % NOTE: the 'actual_split' function uses random init for W and H, so
% %       comparisons will not be exact after the first call to 
% %       actual_split.  Modify 'actual_split' to load a new set of
% %       W and H initializers on each call, as an option.

[tree, splits, is_leaf, clusters, timings, Ws, priorities] = hier8_neat(A, k, params, Winit, Hinit);

