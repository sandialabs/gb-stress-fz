%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Binning as a function of parameter p and determining quartiles of  
% weighted distributions
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% Inputs:
% var           variable to bin and get quartiles of
% p             parameter to perform the binning
% pedges        edges of the bins
% w             list of weights
%
%%% Outputs:
% var_binned    quartiles of the weighted distribution 
%

function var_binned = binquart(var,p,pedges,w)

assert(all(size(var)==size(p)) && size(var,2)==1,'var and p must column vectors of equal size')
assert(size(var,1)==size(w,1),'w must have the same number of rows as var and p')

wi = prod(w,2);             % Weights
var_binned = nan(length(pedges)-1,5);
for i = 1:(length(pedges)-1)                    % Binning according to pedges
    [varsorted,sortI] = sort(var(p>=pedges(i) & p<pedges(i+1)));                % mSmax
    wsorted = wi(sortI);
    wcum = cumsum(wsorted);
    if isempty(wcum)
        var_binned(i,:) = zeros(1,5);
    else
        wtot = wcum(end);
        for j = 1:5                                 % 0:0.25:1 quartiles
            q = -0.25+j*0.25;                       % Quartile
            qI = find(wcum >= wtot*q,1,'first');    % Quartile index
            var_binned(i,j) = varsorted(qI);        % Quartile value
        end
    end
end