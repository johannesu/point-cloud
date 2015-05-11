function compile_script(base_name, sources, extra_args)


base_path = [fileparts(mfilename('fullpath'))];       

if nargin < 2
	sources = {};
end

if nargin < 3
	extra_args = {};
end


if ispc
    % Windows settings.
    
    %%%%%%%%%
    % ACTION REQUIRED:
    % Point this to your Eigen include directory.
    extra_args{end+1} = '-I"C:\Program Files\Eigen"';
    %%%%%%%%

    % Default installation directories for Spii
    extra_args{end+1} = '-I"C:\Program Files\SPII\include"';
    extra_args{end+1} = '-L"C:\Program Files\SPII\lib"';
else
    % Linux settings.
    extra_args{end+1} = '-L/usr/local/lib';
    extra_args{end+1} = '-I/usr/local/include/spii-thirdparty';
    extra_args{end+1} = '-I/usr/local/include/eigen3';
    extra_args{end+1} = '-I/usr/local/include/spii';
    
    extra_args{end+1} = '-lgomp';
end

if ~ispc
    CXXFLAGS = '-fopenmp -std=c++0x';
    extra_args{end+1} = ['CXXFLAGS="\$CXXFLAGS ' CXXFLAGS '"'];
else
	% NOMINMAX fixes windows.h if it gets included.
	extra_args{end+1} = '-DNOMINMAX=1';
    extra_args{end+1} = 'COMPFLAGS="$COMPFLAGS /openmp"';
end

% Link with external libaries
extra_args{end+1} = '-lspii';
    
mex_file_name = [base_path filesep base_name '_mex.' mexext];
cpp_file_name = [base_path filesep base_name '_mex.cpp'];

% Cpp file not found returning
if (isempty(cpp_file_name))
    warning('Mex file not found \n');
    return;
end


mex_file = dir(mex_file_name);
cpp_file = dir(cpp_file_name);

if length(mex_file) == 1
    mex_modified = mex_file.datenum;
else
    mex_modified = 0;
end

cpp_modified   = cpp_file.datenum;

% If modified or not existant compile
compile_file = false;
if ~exist(mex_file_name, 'file')
	compile_file = true;
elseif mex_modified < cpp_modified
	compile_file = true;
end

include_folders = {};
for i = 1 : length(sources);
	include_folders{i} = ['-I' base_path filesep fileparts(sources{i}) filesep];
	sources{i} = [base_path filesep sources{i}];
end


%% Compile
if compile_file
    disp(['Compiling ' cpp_file_name '...']);

    mex_argument = {'mex', ...
        cpp_file_name,  ...
        '-outdir', ...
        base_path, ...
        '-largeArrayDims', ...
        extra_args{:},...
        include_folders{:}, ...
        sources{:}};
    
    % Spaceing needed
    for i = 1:numel(mex_argument)
        mex_argument{i} = [mex_argument{i} ' '];
    end
        
    % Using eval to allow for arguments changeing c and c++ flags
    eval([mex_argument{:}])

    disp('done.');
end

