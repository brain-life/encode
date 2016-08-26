%% Script to compare old vs new results
clear
clc
%diary on

vista_soft_path = '/N/dc2/projects/lifebid/code/vistasoft/';
addpath(genpath(vista_soft_path));

% Define path to the NEW LiFE
new_LiFE_path = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/lifebid_to_share/life/';
addpath(genpath(new_LiFE_path));

dataRootPath = '/N/dc2/projects/lifebid/2t1/HCP/';
dataOutputPath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/Single_TC/105115/';

subject = '105115';

conn = 'NUM01'; % 

lparam_set = {'lmax10','lmax12','lmax2','lmax4','lmax6','lmax8'};
alg_set = {'SD_PROB', 'SD_STREAM'};

ETC_fib = [];

p = 1;

%% TENSOR
alg = 'tensor';
param = '';
alfiles(p).alg = alg;
alfiles(p).param = param;
alfiles(p).connectome = conn;
alfiles(p).subject = subject;

fgFileName    = deblank(ls(fullfile(dataRootPath,subject,'fibers_new', strcat('*b2000*',char(param),'*',char(alg),'*',conn,'*','500000.tck'))));

fg = dtiImportFibersMrtrix(fgFileName);
fg_struct.fg = fg;

% random selection of fibers for ETC
rnd_range = randperm(500000);
rnd_range = rnd_range(1:38468);
ETC_fib = [ETC_fib; fg.fibers(rnd_range)];
alfiles(p).ETC_range = rnd_range;

disp(['Gen results ', alfiles(p).subject, '  connectome' ,alfiles(p).connectome,', ',alfiles(p).alg,' ',alfiles(p).param,'...'])
%Gen_results(dataRootPath, dataOutputPath, subject, alg, param, fg_struct, conn);

clear fg fg_struct

p = p + 1;

%% lmax 2 to lmax 12 for PROB and STREAM (12 connectomes)
for param = lparam_set
    param = char(param);
    for alg = alg_set
        alg = char(alg);
        alfiles(p).alg = alg;
        alfiles(p).param = param;
        alfiles(p).connectome = conn;
        alfiles(p).subject = subject;

        fgFileName    = deblank(ls(fullfile(dataRootPath,subject,'fibers_new', strcat('*b2000*',char(param),'*',char(alg),'*',conn,'*','500000.tck'))));

        fg = dtiImportFibersMrtrix(fgFileName);
        fg_struct.fg = fg;

        % random selection of fibers for ETC
        rnd_range = randperm(500000);
        rnd_range = rnd_range(1:38461);
        ETC_fib = [ETC_fib; fg.fibers(rnd_range)];
        alfiles(p).ETC_range = rnd_range;

        disp(['Gen results ', alfiles(p).subject, '  connectome' ,alfiles(p).connectome,', ',alfiles(p).alg,' ',alfiles(p).param,'...'])
        Gen_results(dataRootPath, dataOutputPath, subject, alg, param, fg_struct, conn);

        clear fg fg_struct

        p = p + 1;

    end
end

%% ETC
% dataOutputPath = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Results/ETC_Dec2015/ETC_full_range_lmax/';
% alg = 'ETC';
% param = '';
% 
% alfiles(p).alg = alg;
% alfiles(p).param = param;
% alfiles(p).connectome = conn;
% alfiles(p).subject = subject;
% 
% fg_struct.fg.fibers = ETC_fib;
% fg_struct.fg.name = 'ETC';
% 
% disp(['Gen results ', alfiles(p).subject, '  connectome' ,alfiles(p).connectome,', ',alfiles(p).alg,' ',alfiles(p).param,'...'])
% Gen_results(dataRootPath, dataOutputPath, subject, alg, param, fg_struct, conn);
% 
% disp('SAVING RESULTS...')
% save(fullfile(dataOutputPath,sprintf('ETC_connectome_info_%s_conn%s.mat',subject,conn)), 'alfiles','-v7.3')

rmpath(genpath(new_LiFE_path));
rmpath(genpath(vista_soft_path));







