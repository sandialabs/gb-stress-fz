%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conglomerate results at the specified symmetry level.
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%% Inputs:
% r         Input vector with results
% id        Id matrix
% level     Level of symmetry to map the results to:
%           0 -> combining all conditions
%           1 -> disorientation
%           2 -> boundary plane
% fun       Function to conglomerate the data:
%           @() -> anonymous function
%           'max' -> max()
%           'min' -> min()
%           'mean' -> weighted mean and standard deviation
%           'stats' -> min, max, weighted mean, std, median, quartiles
% w         Weights vector (only needed for fun='mean' or 'stats')
%
%%% Outputs
% rout      Conglomerated vector
% idout     List of rows to reconstruct the structure at the specified 
%           level
%

function [rout,rows] = symmetryresults(varargin)

    % Input variables and asserts
    [r,id,level,fun] = varargin{1:4};
    if strcmp(fun,'mean') || strcmp(fun,'stats')
        assert(length(varargin)==5,"If fun='mean' or 'stats' then a weights matrix must be provided")
        w = varargin{5};
        assert(all(size(w)==size(id)),'w must have the same size as id')
    end
    assert(size(r,2)==1 && size(r,1)==size(id,1),'r must be an nx1 vector and id an nxN matrix')
    assert(floor(level)==level && (level<size(id,2)),'level must be an integer smaller than the number of columns of id')
    assert(isa(fun,'function_handle') || strcmp(fun,'sum') || strcmp(fun,'min') || strcmp(fun,'max') || ...
           strcmp(fun,'mean') || strcmp(fun,'stats'),'fun must be a function handle or a predetermined string')


    rowi = [];
    
    % Unique ids at the specified level
    idu = unique(id(:,1:level),'rows');
    rows = deal(nan(size(idu,1),1));
    
    % Size of rout
    if isa(fun,'function_handle') || ...
       strcmp(fun,'sum') || strcmp(fun,'min') || strcmp(fun,'max')
        rout = deal(nan(size(idu,1),1));
    elseif strcmp(fun,'mean')
        [rout.mean,rout.std] = deal(nan(size(idu,1),1));
    elseif strcmp(fun,'stats')
        [rout.mean,rout.std,rout.median,rout.min,rout.max] = deal(nan(size(idu,1),1));
        rout.quartile = nan(size(idu,1),5);
    end
    
    for i = 1:size(idu,1)
    
        % Values to consider for each point at the specified level
        I = find(all(id(:,1:level)==idu(i,:),2));
    
        % Data conglomeration
        if strcmp(fun,'sum')
            rout(i) = sum(r(I));
        elseif strcmp(fun,'min')
            [rout(i),Ipoint] = min(r(I));
            rowi = I(Ipoint);
        elseif strcmp(fun,'max')
            [rout(i),Ipoint] = max(r(I));
            rowi = I(Ipoint);
        elseif strcmp(fun,'mean')   
            wi = prod(w(I,:),2);    % Weights
            rout.mean(i) = sum(wi.*r(I))/sum(wi);
            rout.std(i) = sqrt(sum(wi.*(r(I)-rout.mean(i)).^2)/sum(wi));
        elseif strcmp(fun,'stats')   
            wi = prod(w(I,:),2);    % Weights
            rout.min(i) = min(r(I));
            rout.max(i) = max(r(I));
            rout.mean(i) = sum(wi.*r(I))/sum(wi);
            rout.std(i) = sqrt(sum(wi.*(r(I)-rout.mean(i)).^2)/sum(wi));

            [rsorted,sortI] = sort(r(I));       % For median and quartiles:
            wsorted = wi(sortI);
            wcum = cumsum(wsorted);
            wtot = wcum(end);
            for j = 1:5                                 % 0:0.25:1 quartiles
                q = -0.25+j*0.25;                       % Quartile
                qI = find(wcum >= wtot*q,1,'first');    % Quartile index
                rout.quartile(i,j) = rsorted(qI);       % Quartile value
            end
            
            rout.median(i) = rout.quartile(i,3);

        else    % Anonymous function
            rout(i) = fun(r(I));
        end
        
        % rows to reconstruct the test matrix at the specified level
        if ~isempty(rowi)
            rows(i) = rowi;
        else
            rows(i) = I(1);
        end
        
    end


end