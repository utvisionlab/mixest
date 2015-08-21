%% |mxe_readdata|
% *Note:* This is a private function.
%
% Read data points, data weights, indices of data to be worked on, and data
% tag, from the given data argument. (Used in functions that have a data
% input argument)
%
% *Syntax*
%
%   data_struct = mxe_readdata(data)
%   data_struct = mxe_readdata(data, false)
%
% *Description*
%
% The output |data_struct| has the following fields:
%
% * *|data|* (|n-by-N|, where |n| is the data dimensions and |N| is the
% number of total data points) : Data matrix.
% * *|index|* (|1-by-L| where |L| is the size of the demanded subset of
% data points. default value: |[]|) : Data indexing vector for selecting a
% subset of data (Only numeric indexing is supported. Logical indexing is
% not supported.)
% * *|weight|* (|1-by-N|, if index is empty, |1-by-L| if indexing is used.
% default value: |[]|) : Data weight vector,
% * *|size|* (integer): the number of data points (in the demanded subset).
% * *|tag|* (string, default value: |''|): A tag for the data that can be
% used for special purposes.
%
% |data_struct = mxe_readdata(data)| reads the given |data| and fills the
% output structure |data_struct|.
%
% The given |data| may be a matrix, in which case it is returned
% untouched as |data_struct.data|, and the other fields of |data_struct|
% are set to their default values.
%
% In case |data| is a structure, it may contain any of the fields mentioned
% above for |data_struct| (Note that the field |data.data| is mandatory).
% The fields in |data| are copied to the output structure. If |data|
% includes a non-empty |index| field, this syntax applies the indexing to
% the data (|data_struct.data = data.data(:,data.index)|) and empties the
% index (|data_struct.index = []|).
%
% |data_struct = mxe_readdata(data, false)| does not apply the indexing.
% Instead, returns the |index| for customised use. If input |data| is a
% structure without an |index| field (or with an empty value for |index|),
% an index refering to all the data is generated.
%
% *Note:* Use this function only when you use |weight| or |index| in your
% function. The passed in data argument can be passed out to other
% functions intact, when |weight| or |index| are not used in your function.
%
% *Note:* In MATLAB versions R2011b and above, instead of the data matrix,
% you can give a |matlab.io.MatFile| object referring to a MAT-file
% containing a variable named |data| with specifications as pointed out
% above.
%

% Copyright 2015 Reshad Hosseini and Mohamadreza Mash'al
% This file is part of MixEst: visionlab.ut.ac.ir/mixest
%
% Contributors:
%  Reshad Hosseini
%  Mohamadreza Mash'al
%
% Change log: 
%

function data_struct = mxe_readdata(data, applyindex)

    if nargin < 2
        applyindex = true;
    end
    
    data_struct.weight = [];
    data_struct.index = [];
    data_struct.tag = '';
    
    if isstruct(data)
        % given data is a struct
        if isfield(data, 'weight')
            data_struct.weight = data.weight;
        end
        if isfield(data, 'index')
            data_struct.index = data.index;
        end
        if isfield(data, 'tag')
            data_struct.tag = data.tag;
        end
        
        if applyindex && ~isempty(data_struct.index)
            
            % applying the given index
            if isa(data.data, 'matlab.io.MatFile')
                % data.data is a MatFile object => read the variable 'data' from it
                data_struct.data = data.data.data(:, data_struct.index);
            else
                % data.data is a matrix
                data_struct.data = data.data(:, data_struct.index);
            end
            % the index will be emptied below
            
        else
            data_struct.data = data.data;
        end
        
    else
        % given data is a matrix
        data_struct.data = data;
    end
    
    % find N (the size of data)
    if isempty(data_struct.index)
        if isa(data_struct.data, 'matlab.io.MatFile')
            N = size(data_struct.data, 'data', 2); % just a syntax difference for MatFile
        else
            N = size(data_struct.data, 2);
        end
        
        % we don't output an empty index when applyindex is false
        if ~applyindex
            data_struct.index = 1:N;
        end
    else
        N = numel(data_struct.index);
        
        % if index is applied, it is no more valid
        if applyindex
            data_struct.index = [];
        end
    end
    data_struct.size = N;
    
    % find dim
    if isa(data_struct.data, 'matlab.io.MatFile')
        data_struct.dim = size(data_struct.data, 'data', 1); % just a syntax difference for MatFile
    else
        data_struct.dim = size(data_struct.data, 1);
    end
end
