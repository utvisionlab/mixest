function install_mixest

% find required paths
this_file = mfilename('fullpath'); % get this file's path and name (without .m)
toolbox_root_path = fileparts(this_file);

% add MixEst folders to search path 
addpath(toolbox_root_path)
addpath(genpath(fullfile(toolbox_root_path, 'mixest')))
addpath(genpath(fullfile(toolbox_root_path, 'thirdparty')))

% build search database for the documentation
builddocsearchdb(fullfile(toolbox_root_path, 'doc', 'public'))

% compile mex files
mex(fullfile(toolbox_root_path, 'mixest', 'auxillary', 'manopt_manifolds', 'positivedefinite', 'sqrtm_triu_real.c'))
mex(fullfile(toolbox_root_path, 'mixest', 'auxillary', 'manopt_manifolds', 'positivedefinite', 'sqrtm_triu_complex.c'))


