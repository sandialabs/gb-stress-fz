%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shortening a structure (dis, gb or tests)
%
% By Fernando D. León-Cázares
% Sandia National Laboratories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Current options include shortening to to a single point for each
% combination of fundamental zones (disorientation, boundary plane, load
% axis), or to a given set of rows.
%
%%% Inputs:
% s             Structure (dis, gb or tests)
% rows          Rows to preserve:
%               - List of rows
%               - 'uniquefz', the first row for each FZ
%
%%% Outputs:
% s             Shortened structure
% rows          Rows preserved
%

function [s,rows] = sh(s,rows)
    
    sfields = fieldnames(s);

    % Extracting the first element of each fundamental zone
    if strcmp(rows,'uniquefz')      
        sfzfields = sfields(cellfun(@(x) startsWith(x,'fz'),sfields));
        fzs = table;
        for i = 1:length(sfzfields)
            fzs.(sfzfields{i}) = s.(sfzfields{i});
        end
        [~,rows] = unique(fzs,'stable','rows');
    end


    % Shortening the structure s
    for i = 1:length(sfields)
        s.(sfields{i}) = s.(sfields{i})(rows,:);
    end
end