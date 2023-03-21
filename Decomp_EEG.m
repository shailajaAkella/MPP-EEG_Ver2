function [MPP,D]  = Decomp_EEG(X,Clust,M,th,f)
% Function that analyzes/decomposes single-channel EEG traces for a given dictionary
% INPUTS:
% X - structure of bandpassed single trial/ multi - trial EEG traces. Size of structure: number of trials X 1. Field name must be set to Trial. 
% Clust - structure of dictionary atoms. Size of structure: 1 X K  
% th - sparsity constraint from denoising

% if not one, convert to structure
n_tr = size(X,1);
if ~isstruct(X)
    X = squeeze(X);
    n_tr = size(X,1);
    X = cell2struct(X,'Trial',n_tr);
end

% Decomposition per trial
MPP = struct();
for i = 1:n_tr
    [MPP_tr,D] = PhEv_nonovp(X(i).Trial,Clust,M, th,f); % decomposition function
    n_det = size(MPP_tr,2);
    for j = 1:n_det
        MPP_tr(j).PhEv = bsxfun(@rdivide,MPP_tr(j).PhEv,sqrt(sum(MPP_tr(j).PhEv.^2)));  % Unit-norm atoms
    end
    MPP(i).Trials = MPP_tr; 
end

end
