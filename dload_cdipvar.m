function [var_cdip] = dload_cdipvar(cdipid,varname,tlims, tres)
% -------------------------------------------------------------------------
% DLOAD_CDIPVAR  Downloads model data from CDIP THREDDS (historic and rt)
% -------------------------------------------------------------------------
%   Syntax:
%      [var_cdip] = dload_cdipvar(ndbcid,varname,tlims,tres)
%
%   Inputs:
%      CDIPID 
%      VARNAME      options: 'Time' or 'Hs' or 'Dp' or 'Tp' or 'Ta' or 'acmSpeed'
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
%       dload_cdipvar(067,'Hs',[735965 736330],1) 
%       NOTE: this is from 2015 to 2016 every index for the San Nic buoy
%
% Updated as of 08-18-2019 by Alli Ho
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
if varname(1)=='g' | varname(1)=='w' | varname(1:2)=='dw' | varname(1)=='a' | varname(1:2)=='ss'
    paramhead = '';
end

%% LOAD HISTORIC CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
%%% observations (buoys)
mode = 'hist';
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[tstart, tend, starti, endi, ~] = get_tinfo(baseurl, cdipname, tlims,mode, varname);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create URL name and load data
if starti~=endi
    nameurl = [cdipname 'p1/' cdipname 'p1_historic.nc.ascii?'];
    paramurl = [paramhead varname tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    waveVar = dload(url);
    var_cdip_1 = waveVar;
    
    tlims(1) = tend; 
    
else
    var_cdip_1 = [];
    disp('---> Not using CDIP historic');
end

%% LOAD REALTIME CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
%%% observations (buoys)
mode = 'rt';
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[~, ~, starti, endi, ~] = get_tinfo(baseurl, cdipname, tlims,mode,varname);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
% endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% SECOND: create URL name and load data
if starti~=endi & endi>starti;
    nameurl = [cdipname 'p1_rt.nc.ascii?'];
    paramurl = [paramhead varname tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    waveVar = dload(url);
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
if contains(varname, 'Time')
    toff = datenum(1970,1,1,0,0,0);
    var_cdip = var_cdip./(24*60*60) + toff;
end

end

function [tstart tend starti endi Nt] = get_tinfo(baseurl, cdipname, tlims,mode, varname)
    tvarstr = 'waveTime';
    if contains(varname, 'acm')
        tvarstr = 'acmTime';
        % if strcmp(mode, 'hist');
        %     tstart = tlims(1); tend = tlims(2); starti=1; endi=1; Nt=1;
        %     return
        % end
    elseif contains(varname, 'gps')
        tvarstr = 'gpsTime';
    end
    
    if strcmp(mode, 'hist')
        infourl = [baseurl cdipname 'p1/' cdipname 'p1_historic.nc.html'];
    elseif strcmp(mode, 'rt')
        infourl = [baseurl cdipname 'p1_rt.nc.html'];
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
    %%% FIND VAR = Nt   
    idx = find(contains(info, ['Int32 ' tvarstr '[' tvarstr ' = '])); 
    if isempty(idx); tstart = tlims(1); tend = tlims(2); starti=1; endi=1; Nt=1; return; end
    idx = idx(1);    
    
    infot = strsplit(info{idx}, '= '); infot = infot{2};
    Nt = str2num(infot(1:end-2))-1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = tstart, tend
    idx = find(contains(info, 'time_coverage_start'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10));
    tstart = datenum(yr,mo,day);
    
    idx = find(contains(info, 'time_coverage_end'));
    tstr = strsplit(info{idx},' '); 
    tstr = tstr{2};
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10)); hr = str2num(tstr(12:13)); mn = str2num(tstr(15:16));
    tend = datenum(yr,mo,day, hr, mn,0);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     disp(['     CDIP runs from ' datestr(tstart) ' to '  datestr(tend)])
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = starti, endi
    if ~isempty(tlims)
        endi = Nt;
        sdi = 100; % change this if you want to make the time cuttoff more accurate
        starti = endi-floor(endi/sdi)*sdi;
        % downloading coarse resolution time series data to get approximate index location for tlims
        [simpletime, idx] = get_simpletime(baseurl, cdipname, tvarstr, mode, starti, endi, sdi);
        
        [~, idxstart] = min(abs(simpletime-tlims(1)));
        [~, idxend] = min(abs(simpletime-tlims(2)));
        if strcmp(tvarstr, 'acmTime') %idxstart~=idxend
            [ssv, ss] = sort(abs(simpletime-tlims(1)));
            idxstart = min(ss(1:2));
            [ssv, ss] = sort(abs(simpletime-tlims(2)));
            idxend = max(ss(1:2));
        end
        starti = idx(idxstart); endi = idx(idxend);
%         datestr(simpletime([idxstart idxend]))
        if starti~=endi  
            di = 50;
%             di = 20;
            % fine tune starti
            mi = starti;
            sti = mi - di;  if sti<1; sti = 1; end
            eni = mi + di;  if eni>Nt; eni = Nt; end
            sdi = 1;
            [simpletime, idx] = get_simpletime(baseurl, cdipname, varname, mode, sti, eni, sdi);           
            [~, idxstart] = min(abs(simpletime-tlims(1)));
            if strcmp(tvarstr, 'acmTime') %idxstart~=idxend
                [ssv, ss] = sort(abs(simpletime-tlims(1)));
                idxstart = min(ss(1:2));
            end
            starti = idx(idxstart); 
            
            % fine tune endi
            mi = endi;
            sti = mi - di;  if sti<1; sti = 1; end
            eni = mi + di;  if eni>Nt; eni = Nt; end
            sdi = 1;
            [simpletime, idx] = get_simpletime(baseurl, cdipname, varname, mode, sti, eni, sdi);           
            [~, idxend] = min(abs(simpletime-tlims(2)));
            if strcmp(tvarstr, 'acmTime') %idxstart~=idxend
                [ssv, ss] = sort(abs(simpletime-tlims(2)));
                idxend = min(ss(1:2));
            end
            endi = idx(idxend); 
            
        end
        
    else % or small subset (data way to big to pull entire thing)
        starti = Nt-1000;
        endi = Nt;
    end
end

function [simpletime, idx] = get_simpletime(baseurl, cdipname, varname, mode, starti, endi, sdi)
    tvarstr = 'waveTime';
    if contains(varname, 'acm')
        tvarstr = 'acmTime';
    end
    if contains(varname, 'gps')
        tvarstr = 'gpsTime';
    end
    
    %5B43:100:246843%5D
    tstring = [ '[' num2str(starti) ':' num2str(sdi) ':' num2str(endi) ']']; 
    tstring = [ '%5B' num2str(starti) ':' num2str(sdi) ':' num2str(endi) '%5D']; 
    idx = [1 starti:sdi:endi];

    if strcmp(mode, 'hist')
        nameurl = [cdipname 'p1/' cdipname 'p1_historic.nc.ascii?'];
    elseif strcmp(mode, 'rt')
        nameurl = [cdipname 'p1_rt.nc.ascii?'];
    end        
    paramurl = [tvarstr tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    waveVar = dload(url);

    tstring = '[1]';
    tstring = '%5B1%5D';
    paramurl = [tvarstr tstring]; % to do figure out subset!
    url = [baseurl nameurl paramurl];
    initwaveVar = dload(url);

    waveVar = [initwaveVar waveVar];

    toff = datenum(1970,1,1,0,0,0);
    simpletime = double(waveVar)./(24*60*60) + toff; % the course resolution data, corresponding indexes in VAR=idx
    
%     disp(['     simpletime runs from ' datestr(simpletime(1)) ' to '  datestr(simpletime(end))])
end


function [var] = dload(url) % download thredds data from url
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    var = strsplit(data{dstart+1}, ', ');
    var = cellfun(@str2double,var);
end