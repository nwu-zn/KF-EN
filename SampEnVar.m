function [saen,saenvar] = SampEnVar( dim, r, data, tau )
% SAMPEN Sample Entropy
%   calculates the sample entropy of a given time series data

%   SampEn is conceptually similar to approximate entropy (ApEn), but has
%   following differences:
%       1) SampEn does not count self-matching. The possible trouble of
%       having log(0) is avoided by taking logarithm at the latest step.
%       2) SampEn does not depend on the datasize as much as ApEn does. The
%       comparison is shown in the graph that is uploaded.

%   dim     : embedded dimension
%   r       : tolerance (typically 0.2 * std)
%   data    : time-series data
%   tau     : delay time for downsampling (user can omit this, in which case
%             the default value is 1)
%
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Mar 21, 2012
%---------------------------------------------------------------------

if nargin < 4, tau = 1; end
if tau > 1, data = downsample(data, tau); end

N = length(data);
correl = zeros(1,2);
dataMat = zeros(dim+1,N-dim);
for i = 1:dim+1
    dataMat(i,:) = data(i:N-dim+i-1);
end

K=[];   
K_X=[];
for m = dim:dim+1
    count = zeros(1,N-dim);
    K_count = zeros(1,N-dim);
    tempMat = dataMat(1:m,:);
    
    for i = 1:N-m
        % calculate Chebyshev distance, excluding self-matching case
        dist = max(abs(tempMat(:,i+1:N-dim) - repmat(tempMat(:,i),1,N-dim-i)));
        
        ori = tempMat(:,i+1:N-dim);
        obj = repmat(tempMat(:,i),1,N-dim-i);
        second_obj = obj;
        third_obj  = obj;

        if m == 2
            second_obj([1, 2], :) = second_obj([2, 1], :);
            first = abs(ori-obj);
            second = abs(ori-second_obj);
            overlab = min(vertcat(first, second));
        else
            second_obj([1, 2, 3], :) = second_obj([2,3,1], :);
            third_obj([1, 2, 3], :)  = third_obj([3,1,2], :);
            first = abs(ori-obj);
            second = abs(ori-second_obj);
            third = abs(ori-third_obj);
            overlab = min(vertcat(first, second, third));
            
        end
        % calculate Heaviside function of the distance
        % User can change it to any other function
        % for modified sample entropy (mSampEn) calculation
        
        D = (dist < r);
        lab = (overlab < dim);
        lab_count = lab & D;
        
        count(i) = sum(D)/(N-dim);
        K_count(i) = sum(lab_count)/(N-dim);
    end
   
    correl(m-dim+1) = sum(count)/(N-dim);
    K=[K sum(count)];%
    K_X = [K_X sum(K_count)];
end


saen = log(correl(1)/correl(2));

A=K(2);
B=K(1);
KA = K_X(2);
KB = K_X(1);

CP = A/B;
saenvar = abs(CP.*(1-CP)./B)+ 1/(B.^2).*(KA-KB.*CP.^2);
end

