function fig=main(varargin)
folder=fileparts(mfilename('fullpath'));
addpath(folder,'-begin');
fig=app(varargin{:});
end
