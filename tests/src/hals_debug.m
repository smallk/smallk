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
[A] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
               'rnd_256_256.csv']);

[W] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
               'w_init_256_16.csv']);
[H] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
               'h_init_16_256.csv']);

% sparse test
% [A] = mmread(['/Users/richardboyd/Downloads/nmf_package/' ...
%               'reuters.mtx']);
% [W] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
%                'reuters_w_init_128.csv']);
% [H] = csvread(['/Users/richardboyd/Downloads/nmf_package/' ...
%                'reuters_h_init_128.csv']);

fprintf('dimensions of A: %d x %d\n', size(A, 1), size(A, 2));
fprintf('dimensions of W: %d x %d\n', size(W, 1), size(W, 2));
fprintf('dimensions of H: %d x %d\n', size(H, 1), size(H, 2));

fprintf('norm of W: %.15f\n', norm(W, 'fro'));
fprintf('norm of H: %.15f\n', norm(H, 'fro'));

m = size(A, 1);
n = size(A, 2);
k = size(H, 1);

fprintf('m: %d, n: %d, k: %d\n', m, n, k);

left = H * H';
right = A * H';

fprintf('norm of HHt: %.15f\n', norm(left, 'fro'));
fprintf('norm of AHt: %.15f\n', norm(right, 'fro'));

MAXITER = 5000;
tol = 1.0e-4;

for iter = 1 : MAXITER
    
    fprintf('iteration: %d\n', iter);
    
    for i = 1 : k
        W(:, i) = W(:, i) + (right(:, i) - W * left(:, i)) / left(i, i);
        W(:, i) = max(W(:, i), 0);
        if length(find(W(:, i))) == 0
            W(:, i) = eps;
        end
        %fprintf('[%d]\tNorm of W(:, i): %.15f\n', i, norm(W(:,i)));
        W(:, i) = W(:, i) ./ norm(W(:, i));
    end
    
    fprintf('norm of W: %.15f\n', norm(W, 'fro'));
    
    left = W' * W;
    right = A' * W;
    
    fprintf('norm of WtW: %.15f\n', norm(left, 'fro'));
    fprintf('norm of AtW: %.15f\n', norm(right, 'fro'));
    
    for i = 1 : k
        H(i, :) = H(i, :) + (right(:, i)' - left(:, i)' * H) / left(i, i);
        H(i, :) = max(H(i, :), 0);
    end
    
    fprintf('norm of H: %.15f\n', norm(H, 'fro'));
    
    gradH = left * H - right';
    fprintf('norm of gradH: %.15f\n', norm(gradH, 'fro'));
    
    left = H * H';
    fprintf('norm of HtH: %.15f\n', norm(left, 'fro'));
    
    right = A * H';
    fprintf('norm of AHt: %.15f\n', norm(right, 'fro'));
    
    gradW = W * left - right;
    fprintf('norm of gradW: %.15f\n', norm(gradW, 'fro'));

    
	if iter == 1
        initgrad = sqrt(norm(gradW(gradW<=0|W>0))^2 + ...
                        norm(gradH(gradH<=0|H>0))^2);
        fprintf('initgrad value: %.15f\n', initgrad);
		continue;
	else
        projnorm = sqrt(norm(gradW(gradW<=0|W>0))^2 + ...
                        norm(gradH(gradH<=0|H>0))^2);            
     end
 
     progress_metric = projnorm / initgrad;
     fprintf('\tprogress metric: %.15f\n', progress_metric);
 
     if projnorm < tol * initgrad
         fprintf('converged on iteration %d\n', iter);
         break;
     end
end

        
