%% |mxe_htmltable|
% *Note:* This is a private function.
%
% Generates HTML code for a table corresponding to the given structure
% array. (used in documentation)
%
% *Syntax*
%
%   html = mxe_htmltable(S)
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

function html = mxe_htmltable(S)

    html = [...
'<html>' ...
    '<div class="table-responsive">' ...
        '<table class="table table-striped table-bordered">' ...
            '<thead>' ...
                '<tr>' ...
        ];
    
    % table header
    names = fieldnames(S);
    for i = 1 : numel(names)
        html = [html sprintf('<th>%s</th>', strrep(names{i},'_',' ') )]; %#ok<AGROW>
    end
    
    html = [html ...
            '</tr>' ...
        '</thead>' ...
        '<tbody>' ...
        ];

    % table rows
    for n = 1 : numel(S)
        html = [html '<tr>']; %#ok<AGROW>
        for i = 1 : numel(names)
            html = [html sprintf('<td><tt>%s</tt></td>', value2str( S(n).(names{i}) ) )]; %#ok<AGROW>
        end
        html = [html '</tr>']; %#ok<AGROW>
    end
    
    html = [html ...
            '</tbody>' ...
        '</table>' ...
    '</div>' ...
'</html>' ...
        ];
end

function str = value2str(value)

    str = value;
    if isnumeric(value)
        if isempty(value)
            str = '[]';
        else
            str = num2str(value);
        end
    elseif islogical(value)
        if str
            str = 'true';
        else
            str = 'false';
        end
    end
end
