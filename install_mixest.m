function install_mixest

fprintf('Installing the MixEst toolbox...\n');

% Find required paths
this_file = mfilename('fullpath'); % get this file's path and name (without .m)
toolbox_root_path = fileparts(this_file);

% Add MixEst folders to search path
fprintf('Adding MixEst folders to search path...\n');
addpath(toolbox_root_path)
addpath(genpath(fullfile(toolbox_root_path, 'mixest')))
addpath(genpath(fullfile(toolbox_root_path, 'thirdparty')))
fprintf('[Done]\n')

% Compile mex files
fprintf('Compiling mex files...\n');
mex(fullfile(toolbox_root_path, 'mixest', 'auxillary', 'manopt_manifolds', 'positivedefinite', 'sqrtm_triu_real.c'));
mex(fullfile(toolbox_root_path, 'mixest', 'auxillary', 'manopt_manifolds', 'positivedefinite', 'sqrtm_triu_complex.c'));
fprintf('[Done]\n')

fprintf('\nMixEst has been successfully installed.\n');
fprintf('You can use "Set Path" -> "Save" from MATLAB menus to save the added paths permanently.\n\n');

% build search database for the documentation
fprintf('Building search database for the documentation (Optional Step)...\n')
builddocsearchdb(fullfile(toolbox_root_path, 'doc', 'public'))
fprintf('[Done]\n')

fprintf('\nAll Done!\n\n')
