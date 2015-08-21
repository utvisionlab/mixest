%% mxe_publishdocs
% *Note:* This is a private function.
%
% Publish all documentation files into Jekyll source directory
%
% *Syntax*
%
%   mxe_publishdocs
%   mxe_publishdocs online
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


function mxe_publishdocs(arg)

% Note: following fixes should be applied to MATLAB's publish function 
% ([matlabroot]/toolbox/matlab/codetools/publish.m):
%
%
% In the function "getRenderingFigure" change the 'FontSize' from 22 to 32
% (for better readability of rendered equations):
%
%     temptext = text('Parent',tempaxes,'Position',[.5 .5], ...
%         'HorizontalAlignment','center','FontSize',22, ...
%         'Interpreter','latex');
%
%
% in the function "setCodeToEvaluateIfEmpty" remove the following check
% (due to MATLAB 2015a change of "fullfile" behavior with double-dots):
%
%     if ~strcmpi(strrep(fullPathToScript,'/',filesep),which(cmd)) && ...
%         (options.evalCode==true)
%         error(pm('OffPath'))
%     end
%
%


    if nargin>0 && strcmpi(arg, 'online')
        online = true;
    else
        online = false;
    end
    
    % find required paths
    this_file = mfilename('fullpath'); % get this file's path and name (without .m)
    this_file_path = fileparts(this_file);
    
%     toolbox_root_path = fullfile(this_file_path, '..', '..', '..');
    % (MATLAB 2015a changed behavior of "fullfile" makes "publish" unhappy)
    toolbox_root_path = this_file_path;
    for k = 1:3
        toolbox_root_path = fileparts(toolbox_root_path);
    end
    
    doc_source_path = fullfile(toolbox_root_path, 'doc', 'source');
    xsl_path = fullfile(doc_source_path, '_xsl');
    if online
        publish_dest_path = fullfile(toolbox_root_path, '..', 'jekyll_online', 'docs');
        main_rootpath = '../';
    else
        publish_dest_path = fullfile(doc_source_path, '_jekyll');
        main_rootpath = '';
    end

    % build common publish options
    publish_options.createThumbnail = false;
    publish_options.stylesheet = fullfile(xsl_path, 'mxdom2jekyll.xsl');
    
    % publish main documentation files
    publish_options.evalCode = true;
    publish_options.showCode = false;
    fprintf('\nPublishing main documentation files...\n')
    publishall(doc_source_path, publish_dest_path, true, publish_options, main_rootpath, 1, '*');

    % publish API reference
    publish_options.evalCode = false;
    publish_options.showCode = false;
    api_source_path = fullfile(toolbox_root_path, 'mixest');
    api_dest_path = fullfile(publish_dest_path, 'api');
    api_rootpath = [main_rootpath '../'];
    fprintf('\nPublishing API reference...\n')
    index_html = publishall(api_source_path, api_dest_path, false, publish_options, api_rootpath, 1, '*.m');
    % add an index.html for the API reference
    jekyll_front_matter = sprintf([ ...
        '---\n' ...
        'rootpath: "%s"\n' ...
        'layout: page\n' ...
        'title: "%s"\n' ...
        '---\n' ...
        ], ...
        api_rootpath, 'API Reference');
    index_html = [jekyll_front_matter index_html];
    dlmwrite(fullfile(api_dest_path, 'index.html'), index_html, 'delimiter', '');
    
    % publish examples
    publish_options.evalCode = false;
    publish_options.showCode = true;
    examples_source_path = fullfile(toolbox_root_path, 'examples');
    examples_dest_path = fullfile(publish_dest_path, 'examples');
    examples_rootpath = [main_rootpath '../'];
    fprintf('\nPublishing examples...\n')
    index_html = publishall(examples_source_path, examples_dest_path, false, publish_options, examples_rootpath, 1, 'example*.m');
    % add an index.html for the examples
    jekyll_front_matter = sprintf([ ...
        '---\n' ...
        'rootpath: "%s"\n' ...
        'layout: page\n' ...
        'title: "%s"\n' ...
        '---\n' ...
        ], ...
        examples_rootpath, 'Examples');
    index_html = [jekyll_front_matter index_html];
    dlmwrite(fullfile(examples_dest_path, 'index.html'), index_html, 'delimiter', '');

    fprintf('\nAll Done!\n\n')
end

function index_html = publishall(src_path, dest_path, replicate_dirs, publish_options, rootpath, depth, wildcard)
% recursive function to publish all M-files from the given source path 
% and all its sub-directories to the given destination path.
%
% src_path: source doc path
% dest_path: publish output path
% replicate_dirs: when true, the directory structure of the source is
%   replicated at the destination
% publish_options: publish options except outputDir which will be
%   overwritten
% rootpath: relative path to the root directory where the css, js, ...
%   directories required by the produced html reside.
% depth: current depth in folders
% index_html = publishall(...): returns a string containing links to 
%   published output files as html tags.

    publish_options.outputDir = dest_path;
    index_html = '';
    
    % first publish the files in src_path
    listing = dir(fullfile(src_path, wildcard));
    dirIdx = [listing.isdir];
    file_listing = listing(~dirIdx);
    excludeIdx = ismember({file_listing.name}, {'Contents.m', 'info.xml'}); % exclude special files
    if any(excludeIdx)
        file_listing = file_listing(~excludeIdx);
    end
    if numel(file_listing) > 0
        if publish_options.evalCode
            oldwd = cd(src_path); % change directory to make publish happy
        end
        index_html = '<ul>';
        for k = 1:numel(file_listing)
            file_name = file_listing(k).name;
            if ~ismember(file_name(1), {'.', '_'}) % don't process any file name starting with '.' or '_'
                [unused1, unused2, ext]= fileparts(file_name); %#ok<ASGLU>
                if strcmpi(ext, '.m')
                    % m-files are published
                    src_file = fullfile(src_path, file_name);
                    fprintf('Publishing: %s ', src_file)
                    published_doc = publish(src_file, publish_options);
                    addFrontMatter(published_doc, rootpath);
                    fprintf('-> %s [Done]\n', published_doc)

                    [unused, func_name]= fileparts(file_name); %#ok<ASGLU>
                    index_html = [index_html ...
                        sprintf('<li><a href="%s.html">%s</a></li>', func_name, func_name), ...
                        ]; %#ok<AGROW>
                else
                    % other files are copied directly
                    src_file = fullfile(src_path, file_name);
                    dest_file = fullfile(dest_path, file_name);
                    fprintf('   Copying: %s -> %s ', src_file, dest_file)
                    copyfile(src_file, dest_file);
                    if strcmpi(ext, '.html')
                        % html files' rootpath should also be updated
                        addFrontMatter(dest_file, rootpath);
                    end
                    fprintf('[Done]\n')
                end
            end
        end
        index_html = [index_html '</ul>'];
        if publish_options.evalCode
            cd(oldwd);
        end
    end

    % now, check for any sub-directories
    listing = dir(src_path);
    dirIdx = [listing.isdir];
    dir_listing = listing(dirIdx);
    for k = 1:numel(dir_listing)
        dir_name = dir_listing(k).name;
        if ~ismember(dir_name(1), {'.', '_'}) % don't process '.', '..' or any directory name starting with '.' or '_'
            if replicate_dirs
                new_dest_path = fullfile(dest_path, dir_name);
                new_rootpath = [rootpath '../'];
                if ~exist(new_dest_path, 'dir')
                    mkdir(new_dest_path)
                end
            else
                new_dest_path = dest_path;
                new_rootpath = rootpath;
            end
            s = publishall(fullfile(src_path, dir_name), new_dest_path, ...
                replicate_dirs, publish_options, new_rootpath, depth+1, wildcard);
            index_html = [index_html ...
                sprintf('<h%d>%s</h%d>', depth+1, dir_name, depth+1) ...
                s]; %#ok<AGROW>
        end
    end
    
end

function addFrontMatter(published_doc, rootpath)
% add the first lines of Jekyll front matter to the published html
% document. Note: the last lines containing the document title are defined
% in the publishing XSL file (mxdom2jekyll.xsl).

    frontMatter = sprintf([ ...
        '---\n' ...
        'rootpath: "%s"\n' ...
        ], ...
        rootpath);

    % this is how we should insert lines at the beginning of a file in MATLAB
    docText = fileread(published_doc);
    dlmwrite(published_doc, [frontMatter docText], 'delimiter', '');
end