function [var_cdip, t_cdip] = dload_cdipxyz(cdipid,varname,tlims, tres)
% -------------------------------------------------------------------------
% DLOAD_CDIPXYZ  Downloads CDIP displacement data from CDIP THREDDS (historic and rt)
% -------------------------------------------------------------------------
%   Syntax:
%      [var_cdip] = dload_cdipxyz(cdipid,varname,tlims,tres)
%
%   Inputs:
%      CDIPID 
%      VARNAME      options: 'XDisplacement' or 'YDisplacement' or 'ZDisplacement'
%      TLIMS        default = [Nt-100 to Nt];
%      TRES         default = 1;
% 
%   Output:
%      [var_cdip]
% 
%   Uses:
%      get_tinfo.m (subfunction)
% 
%   Sample: 
% 
%   Other info:
%       see http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/220p1_xy.nc.html
%       for varnames
%
% Updated as of 03-15-2023 by Alli Ho
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

%% function set up 
if isnumeric(cdipid)
    cdipid = num2str(cdipid);
    if length(cdipid)<3
        cdipid = ['0' cdipid];
    end
end
cdipname = cdipid;

if ~exist('tlims')
    tlims = [];
end

if ~exist('tres')
    tres = 1;
end

paramhead = 'wave';
if varname(1)=='g' | varname(1)=='w' | varname(1:2)=='dw' | varname(1)=='a'
    paramhead = '';
end
paramhead = 'xyz';
var_cdip = []; t_cdip = [];
%% LOAD HISTORIC RAW DEPLOYMENT DATA CDIP DATA FROM CDIP THREDDS
% % -------------------------------------------------------------------------
%%% observations (buoys)
mode = 'dep';
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get list of deployments
[depnums] = get_allxyzdeployments(baseurl, cdipname);
[depnums] = choose_xyzdeployments(baseurl, cdipname, tlims, depnums);
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% get data from each deployment

var_cdip_1 = [];
t_cdip_1 = [];
if ~isempty(depnums)
    for i=1:length(depnums)
        depnum = depnums(i);
        depstr = num2str(depnum);
        if length(depstr)<2
            depstr = ['0' depstr];
        end
        %%% FIRST: get tinfo and limit time series if need be
        [tstart, tend, starti, endi, ~, dt] = get_tinfo(baseurl, cdipname, tlims,mode, depnum);
        % check if starti, endi, tres, are even/logical
        tt = [starti:tres:endi];
        endi = tt(end);
        tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
        tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 

        %%% SECOND: create URL name and load data
        if starti~=endi
            nameurl = [cdipname 'p1/' cdipname 'p1_d' depstr '.nc.ascii?'];
            paramurl = [paramhead varname tstring]; % to do figure out subset!
            url = [baseurl nameurl paramurl];
            waveVar = dload(url);
            var_cdip_1 = [var_cdip_1 waveVar];

            tt = [1:tres:endi-starti+1];
            waveVar = tlims(1):dt:tlims(2)+dt;
            waveVar = waveVar(tt);
            t_cdip_1  = [t_cdip_1 waveVar];

            tlims(1) = tend; 

        end
    end
else
    disp('---> Not using CDIP historic deployment data');
end

% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

var_cdip_2 = [];
t_cdip_2 = [];
var_cdip = [var_cdip_1 var_cdip_2];
t_cdip = [t_cdip_1 t_cdip_2];


%% LOAD REALTIME CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
%%% observations (buoys)
mode = 'rt';
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/';
%baseurl = 'https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/132p1_xy.nc.html
% .../132p1_xy.nc.ascii?xyzXDisplacement[0:1:200]
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[tstart, tend, starti, endi, ~, dt] = get_tinfo(baseurl, cdipname, tlims,mode);
if tstart>tlims(1)
%     var_cdip = []; t_cdip = [];
    disp('---> Not using CDIP realtime');
    return
end
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
% endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create URL name and load data
if starti~=endi & endi>starti;
    nameurl = [cdipname 'p1_xy.nc.ascii?'];
    paramurl = [paramhead varname tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    waveVar = dload(url);
    waveVar(waveVar==-999.99) = NaN;
    var_cdip_2 = waveVar;
else
    var_cdip_2 = [];
    disp('---> Not using CDIP realtime');
end

%% COMBINE REALTIME AND HISTORIC
% -------------------------------------------------------------------------
var_cdip = [var_cdip_1 var_cdip_2];

%% IF TIME TURN TO DATENUM FORMAT
% -------------------------------------------------------------------------
tt = [1:tres:endi-starti+1];
t_cdip = tlims(1):dt:tlims(2)+dt;
t_cdip = t_cdip(tt);

end

function [depnums] = get_allxyzdeployments(baseurl, cdipname)
% xyz displacements and waveCheckFactor not included

    listurl = ['http://thredds.cdip.ucsd.edu/thredds/catalog/cdip/archive/' cdipname 'p1/catalog.html'];
    listinfo = urlread(listurl); listinfo = strsplit(listinfo,'\n');
    idx = find(contains(listinfo, ['p1_d']));
    depnums = cellfun(@(x) strsplit(x, 'tt>'), listinfo(idx), 'UniformOutput', false);
    depnums = cellfun(@(x) x(8:end-5), cellfun(@(x) x{2}, depnums, 'UniformOutput', false), 'UniformOutput', false);
    
    DEPIDX = ones(size(depnums));
    for i=1:length(depnums)
        depnum = str2num(depnums{i});
        depstr = num2str(depnum);
        if length(depstr)<2
            depstr = ['0' depstr];
        end
        infourl = [baseurl cdipname 'p1/' cdipname 'p1_d' depstr '.nc.html'];
        info = urlread(infourl); info = strsplit(info,'\n');
        idx = find(contains(info, ['xyz displacements and waveCheckFactor not included'])); 
        if ~isempty(idx); DEPIDX(i) = 0; end
    end
    ii = find(DEPIDX);
    depnums = cellfun(@(x) str2num(x), depnums(ii)); 
end

function [depnums] = choose_xyzdeployments(baseurl, cdipname, tlims, depnums)
    DEPTLIMS = NaN(2,length(depnums));
    for i=1:length(depnums)
        depnum = depnums(i);
%         depstr = num2str(depnum);
%         if length(depstr)<2
%             depstr = ['0' depstr];
%         end
%         infourl = [baseurl cdipname 'p1/' cdipname 'p1_d' depstr '.nc.html'];
%         info = urlread(infourl); info = strsplit(info,'\n');
        [tstart tend starti endi Nt dt] = get_tinfo(baseurl, cdipname, tlims,'dep',depnum);
        DEPTLIMS(:,i) = [tstart tend];
        
    end
    iis = find(DEPTLIMS(1,:)<tlims(2)); 
    iie = find(DEPTLIMS(2,:)>tlims(1)); 
    ii = intersect(iie, iis);
    depnums = depnums(ii);
%     iie = find(DEPTLIMS(2,:)<tlims(1)); 
    
end


function [tstart tend starti endi Nt dt] = get_tinfo(baseurl, cdipname, tlims,mode,varargin)
%     https://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/ 132p1_xy.nc.ascii?xyzStartTime,xyzSampleRate
    if strcmp(mode, 'hist')
        infourl = [baseurl cdipname 'p1/' cdipname 'p1_historic.nc.html'];
    elseif strcmp(mode, 'rt')
%         infourl = [baseurl cdipname 'p1_xy.nc.ascii?xyzStartTime,xyzSampleRate'];
        infourl = [baseurl cdipname 'p1_xy.nc.html'];
    elseif strcmp(mode, 'dep')
        depnum = varargin{1};
        depstr = num2str(depnum);
        if length(depstr)<2
            depstr = ['0' depstr];
        end
        infourl = [baseurl cdipname 'p1/' cdipname 'p1_d' depstr '.nc.html'];        
    else
        error('Mode for CDIP data not compatible');
    end
    try
        info = urlread(infourl); info = strsplit(info,'\n');
    catch
        tstart = tlims(1); tend = tlims(2); starti=1; endi=1; Nt=1;
        return
    end

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = tstart, tend
    idx = find(contains(info, ['time_coverage_start'])); idx = idx(1);
    infot = strsplit(info{idx}, ': '); tstr = infot{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10)); hr = str2num(tstr(12:13)); min = str2num(tstr(15:16));
    tstart = datenum(yr,mo,day, hr, min,0);

    idx = find(contains(info, ['time_coverage_end'])); idx = idx(1);
    infot = strsplit(info{idx}, ': '); tstr = infot{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10)); hr = str2num(tstr(12:13)); min = str2num(tstr(15:16));
    tend = datenum(yr,mo,day, hr, min,0);

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = Nt   
    idx = find(contains(info, ['[xyzCount ='])); idx = idx(1);
    infot = strsplit(info{idx}, '..'); infot = infot{2};
    Nt = str2num(infot(1:end-1));
    dt = (tend-tstart)./Nt;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     disp(['     CDIP runs from ' datestr(tstart) ' to '  datestr(tend)])
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = starti, endi
    if ~isempty(tlims)
        
        starti = round((tlims(1)-tstart)./dt);
        endi = Nt-round((tend-tlims(2))./dt);
            
        if endi>Nt;
            endi= Nt;
        end
    else % or small subset (data way to big to pull entire thing)
        starti = Nt-1000;
        endi = Nt;
    end
end

function [var] = dload(url) % download thredds data from url
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    var = strsplit(data{dstart+1}, ', ');
    var = cellfun(@str2double,var);
end