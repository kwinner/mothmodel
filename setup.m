function setup()

root = fileparts(which(mfilename()));
addpath(genpath(root));

curdir = pwd();

disp 'Building discrete-ars-0.1'
cd(sprintf('%s/lib/discrete-ars-0.1', root));
make;
cd(curdir);

disp 'Building dbinom'
cd(sprintf('%s/util', root));
mex dbinom.c
cd(curdir);


end