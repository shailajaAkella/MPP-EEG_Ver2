function [MPP,th] = Denoise(X,M,mul)
% Function to isolate putative neuromodulatory events from background noise
% Uses a hierarchical detection scheme and a correntropic similarity
% measure to detect PUTATIVE events
% INPUTS: 
% X - number of trials X 1 cells/struct comprising of trial signals
% M - maximum length of neuromodulation
% mul - multiplicative factor for Silverman's rule
% USES CORRENTROPY AND STANDARD DEVIATION 
% OUTPUTS: 
% MPP - structure with putative neuromodulatory detections and features
% th - minimum norm below which every M - length snippet is classified as
% background noise

n_tr = size(X,1);

% check if X is a structure and convert if not
if ~isstruct(X)
    X = squeeze(X);
    n_tr = size(X,1);
    X = cell2struct(X,'Trial',n_tr);
end

% initializations
f = []; X_M = []; idtr = [];

% Extracting all M - length snippets from all trials
for i = 1:n_tr
    X_temp1 = X(i).Trial;
    X_temp2 = X_temp1;
    
    env = abs(hilbert(X_temp1));
    aux_M = round(M/2);
    spn = round(aux_M/2)*2 - 1;
    sm = smooth(env,spn);
    
    % First, find clear peaks
    [~,pk_loc] = findpeaks(sm,'MinPeakDistance',M, 'SortStr','descend');
    N = length(X_temp1);
   
    for j = 1:length(pk_loc)
        if pk_loc(j) - round(M/2) <= 0
            t = 1:pk_loc(j) + round(M/2) - 1;
            if length(t) > round(0.75*M)
                X_M = [X_M; [X_temp1(t) X_temp1(t(end) + 1:t(end) + M - length(t))]];
                X_temp1(t) = 0;
                f = [f; t(1)];
                idtr = [idtr i];
            end
        elseif pk_loc(j) + round(M/2) - 1 > N
            t = pk_loc(j) - round(M/2):N;
            if length(t) > round(0.75*M)
                X_M = [X_M; [X_temp1(t(1) - (M - length(t)):t(1) - 1) X_temp1(t)]];
                X_temp1(t) = 0;
                f = [f; t(1)];
                idtr = [idtr i];
            end
        else
            t = round(pk_loc(j) - M/2:pk_loc(j) + M/2 - 1);
            X_M = [X_M; X_temp1(t)];
            X_temp1(t) = 0;
            f = [f; t(1)];
            idtr = [idtr i];
        end
    end
    
    pts = find(X_temp1);
    pt1 = [pts(1) pts(find(diff(pts)>1)+1)];
    pt2 = [pts(find(diff(pts)>1)) pts(end)];
    
    % Second, find remaining M - length snippets
    for k = 1:length(pt2)
        snippet = X_temp2(pt1(k):pt2(k));
        L = length(snippet);
        x = [];
        if (L >= round(0.75*M) && L <= M && pt2(k) + (M-L) <= N)
            x = [snippet X_temp2(pt2(k)+1:pt2(k)+ (M-L))];
            X_M = [X_M; x];
            f = [f; pt1(k)];
            idtr = [idtr i];
        elseif (L > M)
            m = floor(L/M);
            for j = 1:m
                x = [x; snippet(1 +(j-1)*M:M+ (j-1)*M)];
                f = [f; pt1(k)+(j-1)*M];
                idtr = [idtr i];
            end
            if (pt2(k)+(m+1)*M-L <= N && (L - m*M) >= 0.75*M)
                x = [x; X_temp2(pt1(k)+ m*M:pt2(k)+(m+1)*M-L)];
                f = [f; pt1(k) + m*M];
                idtr = [idtr i];
            end
            X_M = [X_M; x];
        end
    end
end

X_M = X_M';
[~,cols] = size(X_M);

CorrMat = zeros(cols,cols);
sig_e = std(X_M(:));
sig = (1.06 * sig_e * (length(X_M(:)))^(-0.2))*mul;

% Third, parallel computation of correntropy matrix
parfor i = 1:cols
    for j = 1:cols
        if (j <= i)
            continue;
        end      
        CorrMat(i,j) =(1/M) * sum((1/(sqrt(2*pi)*sig))*exp((-(X_M(:,i) - X_M(:,j)).^2)/(2*sig^2)));
    end
end

CorrMat = CorrMat + CorrMat';
% similarity vector
P_vec = sum(CorrMat);
P_vec = (1/(cols-1))*(P_vec);
[P_vec,I] = sort(P_vec,'descend'); 

% Four, find putative events - snippets that are one standard deviation below the mean
% P vector value
[~, index] = min(abs(P_vec - (mean(P_vec) - std(P_vec))));

prc = 1 - index/cols; 
prc = round(prc*cols);
idx_sp = I(end - prc:end);

% Last, construct the MPP
MPP = struct();
if (~isempty(idx_sp))
    for i = 1:length(idx_sp)
        MPP(i).PhEv = X_M(:,idx_sp(i));
        MPP(i).norm = norm(X_M(:,idx_sp(i))); 
        MPP(i).tau = f(idx_sp(i)) + round(M/2); 
        MPP(i).idtr = idtr(idx_sp(i));
    end
    th = min([MPP.norm]); % threshold
else
    MPP(1).PhEv = [];
    MPP(1).alpha = [];
    MPP(1).tau = [];
    MPP(1).idtr = [];
end

end


