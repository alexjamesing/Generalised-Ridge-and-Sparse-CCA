function mat_deflate = matrix_deflation(mat,w)
%% function SCCA_matrix_deflation
% The aim of this script is to orthogonalise the data matrix 'mat' with
% respect to previous canonical correlates. 
%
%% Inputs
% mat: This is the data matrix to orthogonalise. We orthogonalise this
% matrix with respect to previous canonical variates. Here, we use the
% canonical weight vector w.
% w: This is the canonical weight vector associated with   

%% Outputs
%  mat_deflate: This is the input matrix, which is now orthogonalised with
%  respect to the canonical variate associated with the canonical weight
%  vector w.
%
%%

mat = zscore(mat);

c = zscore(mat*w);
beta = ((inv(c'*c))*c')*mat;
mat_deflate = zscore(mat-((c*beta)));
    
end
