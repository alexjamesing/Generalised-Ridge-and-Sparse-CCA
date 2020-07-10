
function [c_arr,w_arr,Corr_mat] = GRCCA(X_arr, C_mv,rp_arr,niter)

%% Inputs:
% X_arr: this is a 1 x mv cell array of matrices, where mv is the number of
% matrices to be compared in the analysis conducted. Each cell of the
% matrix should contain a matrix represeting a view of the data.
% rp_arr: this is a 1 x mv array. This array should contain ridge
% values for each of the views to be compared in the analysis.
% w_sgn_arr: this is a 1 x mv array. Value in an array specify whether weight 
% vectors should be constrained to be positive (1), negative(-1) or
% unconstrained (0).
% C_mv: this is a mv x mv matrix. Each element of the matrix specifies
% whether or not the script is designed to optimise relations between those
% data sets.
% niter: number of iterations script must run through
%% Outputs:
% c_arr: this is a 1 x mv cell array of canonical correlates, one for each
% view of the data.
% w_arr: this is a 1 x mv cell array of canonical vectors, one for each
% view of the data.
% Corr_mat: this is a mv x mv correlation matrix. Each element of the matrix
% specifies the correlation between associated canonical correlates.
% 


w_arr = cell(1,size(X_arr,2)); %% this is an array for 'outer' PLS/canoncor coefficients 
COV_arr = cell(1,size(X_arr,2)); %% this is an array of covariance matrices, which can be used later on to speed up calculation
c_arr = zeros(size(X_arr{1,1},1),size(X_arr,2)); %% this is the array in which we store canonical correlate vectors


%% calculate the covariance matrices, for use in later calculations

for i = 1:size(X_arr,2);
    
  Xtemp = X_arr{1,i};   
  dimX = size(Xtemp,2);
  wtemp = randn(dimX,1);
  rp = rp_arr(i);
  COV_arr{i} = inv(rp*(eye(size(Xtemp,2),size(Xtemp,2))) + (1-rp)*cov(Xtemp)); %% calculate the covariance matrices
  w_arr{i} = ((wtemp'*(COV_arr{i})*wtemp)^-0.5)*(COV_arr{i})*wtemp;
 
end



%% compute inner components

  for iter = 1:niter
    iter
    for j = 1:size(X_arr,2); %% in this loop, we apply the msCCA algorithm in each view of the data
          
     w = cca_vec(X_arr,COV_arr,w_arr,C_mv,j); %% this sub-routine actually computes the weights
     w_arr{j} = w;
     c_arr(:,j) = (zscore((X_arr{j})))*w_arr{j}; %% calculate canonical correlates  
     
    end

  end
  
  
Corr_mat = corr(c_arr,c_arr);

end


function w = cca_vec(X_arr,COV_arr, w_arr,C_mv,j) %% this sub-routine is used to calculate the correlation between sets 
 
  w_adder = find(C_mv(j,:));  
 
  
  ca = zeros(size(X_arr{1,1},1),1);
  for x = 1:size(w_adder,2);
  
  X = ((X_arr{w_adder(x)}));    
  ca = ca+(X_arr{w_adder(x)})*(w_arr{w_adder(x)});
    
  end

  X = (X_arr{j});
  ca = ca*((((X'*ca)')*(COV_arr{j})*(X'*ca))^-0.5);
  w = (COV_arr{j})*(X'*ca);    %% compute outer weight coefficient

  
end










