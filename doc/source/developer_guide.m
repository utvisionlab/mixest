%% Developer Guide
% *Note:* This document is a draft.
%
% Here are some hints and guidelines for MixEst developers.
%

%% Coding
% * Take a look at
% <https://www.math.tu-berlin.de/fileadmin/i26/download/AG_ModNumDiff/FG_NumMath/seminars/toolseminar/poloni_slides.pdf
% "Unit tests in Matlab and other good coding practices"> by Federico
% Poloni.
% * Please Don't use the tilde (|~|) operator for unused function outputs
% (as in |[~, ...] = function(...);|). This is not compatible with MATLAB
% versions prior to R2009b. Instead, use the context menu in the MATLAB
% editor to suppress the warning on the line: |[unused, ...] =
% function(...); %#ok<ASGLU>|.
% * Use <api/mxe_readdata.html |mxe_readdata|> before accessing given
% |data|.

%% Estimation functions
% * Always check |options.verbosity| before outputting informative text.
% * Use |mxe_inneroptions| to get options for your inner estimations (if
% any).

%% Structures in place of objects
% * Since we dislike the object-orientation of MATLAB in many ways, we use
% structures created by factory functions in place of objects. This way,
% the object's properties are the variables defined in the workspace of the
% factory function and the object's methods are exposed as fields containig
% a handle to a nested function within the factory function.
% * To give access to read-only properties, use functions (getter methods).
% * Methods that need to change some properties from the current object,
% and return the result as a new object, should call the factory function
% using conventional special inputs. This is required in order to have a
% new workspace for the new properties. See the codes of the |gammafactory|
% and its |fixate| function for an example of using this technique.

%% Documenting
% * We use the |publish| function to generate MixEst documentations. See
% <http://mathworks.com/help/matlab/matlab_prog/marking-up-matlab-comments-for-publishing.html
% Publishing Markup> in MATLAB's documentation to get familiar with the
% markup.
% * A useful shortcut key is |Ctrl+J| which wraps the doc text
% automatically.
% * Please write the documentation for your codes right after coding. The
% descriptions need not be lengthy. Just copy the documentation from a
% similar function and paste it above yours and then change it as needed.
% But never let an irrelevant documentation remain there.
% 

%% Technical details
% * Compiling the documentation is done in two steps: First, by running
% |mxe_publishdocs|, the M-files are published using a special style-sheet
% to prepare them for <http://jekyllrb.com Jekyll>. The published files
% are put in the Jekyll source folder (|mixest/doc/source/_jekyll|). The
% next step is to build them using Jekyll (|jekyll build|).
% * The main documentation files are at |mixest/doc/source| and are
% published to the Jekyll source folder. The API reference is extracted
% from all the M-files under |mixest/mixest| and is all put into a
% sub-directory named |api| beside the published main docs. Therefore, to
% give a link to a page in the API reference from a main documentation, use
% the |api/| prefix in the link like |api/mixturefactory.html|. To give a
% link to a page in the main documentation from a function documentation,
% use the |../| prefix in the link like |../estimation_options.html|.
% * We don't run the codes when publishing. You should add any results
% (like plot images) manually. Put the images in a sub-directory |img|
% beside the source documentation files and add the |img/| prefix to the
% file name in the markup for the image.
%
