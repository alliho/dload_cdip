function [spec_cdip f_cdip t_cdip bw_cdip] = dload_cdipspec(cdipid,varname,tlims,tres)
% -------------------------------------------------------------------------
% DLOAD_CDIPSPEC  Downloads spectral model data from CDIP THREDDS 
% -------------------------------------------------------------------------
%   Syntax:
%      [spec_cdip f_cdip t_cdip] = dload_cdipspec(cdipid,varname,tlims,tres)
%
%   Inputs:
%      CDIPID 
%      VARNAME      options: 'EnergyDensity' or 
%                            'Mean Direction'      
%      TLIMS        default = [NaN NaN], can also use numeric indexes...
%      TRES         default = 1;
% 
%   Output:
%      [spec_cdip f_cdip t_cdip]
% 
%   Uses:
%      get_tinfo.m (subfunction)
% 
% Created 08-18-2019 by Alli Ho
% Updated as of 08-18-2019 by Alli Ho
% Updated 12-13-2021 to account for CDIP server syntax change? 
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
elseif ~isnan(tlims)
    if tlims(2)<tlims(1)
        error('Incorrect time limits')
    end
end

if ~exist('tres')
    tres = 1;
end

if strcmp(varname, 'EnergyDensity') | strcmp(varname, 'MeanDirection') | contains(varname, 'Value') | strcmp(varname, 'CheckFactor') | contains(varname, 'dwr')
    dim = 2;
else
    error('Provided variable name is not real');
end

vo = 0;
fstr = '[0:1:63]'; fN = 64;
fstr = '%5B0:1:63%5D'; fN = 64;




%% LOAD HISTORIC CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
%%% observations (buoys)

mode = 'hist';
% fstr = '%5B0:1:63%5D'; fN = 64;
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/archive/';
nameurl = [cdipname 'p1/' cdipname 'p1_historic.nc.ascii?'];


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[tstart, tend, starti, endi, Nt] = get_tinfo(baseurl, cdipname, tlims,mode);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 

Nt = Nt/tres+1;
if endi~=starti
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% SECOND: get frequency

    % Load frequency
    paramurl = ['waveFrequency' fstr];
    url = [baseurl nameurl paramurl]; %full url name
    f_cdip = dload(url);
    
    % Load bandwidthh
    paramurl = ['waveBandwidth' fstr];
    url = [baseurl nameurl paramurl]; %full url name
    bw_cdip = dload(url);

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% AND THEN: create empty variable, url name, load data, reshape

    %%%% waveEnergyDensity %%%
    if dim==2 
        paramurl = ['wave' varname tstring fstr];
        var_cdip = NaN(fN,Nt);
        url = [baseurl nameurl paramurl]; %full url name
        data = urlread(url); data = strsplit(data,'\n');
        dstart = find(strcmp(data,'---------------------------------------------'))+1;

        for j=1:Nt
            waveVar = strsplit(data{dstart+j}, ', ');
            if vo; disp(waveVar{1}); end
            waveVar = cellfun(@str2double,waveVar);
            var_cdip(:,j) = waveVar(2:end);
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    spec_cdip_1 = var_cdip;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% THIRD: provide time variable
    tstart = find(strcmp(data,['wave' varname '.waveTime[' num2str(Nt) ']']));
    waveVar = strsplit(data{tstart+1}, ', ');
    waveVar = cellfun(@str2double,waveVar);
    toff = datenum(1970,1,1,0,0,0);
    t_cdip_1 = waveVar./(24*60*60) + toff;
    
    tlims(1) = tend;
else
    spec_cdip_1 = [];
    t_cdip_1 = [];
    disp('---> Not using CDIP historic');
end

%% LOAD REALTIME CDIP DATA FROM CDIP THREDDS
% -------------------------------------------------------------------------
if contains(cdipid, '238')
    fstr = '%5B0:1:99%5D'; fN = 100;
end

%%% observations (buoys)
mode = 'rt';
% fstr = '%5B0:1:99%5D'; fN = 100;
baseurl = 'http://thredds.cdip.ucsd.edu/thredds/dodsC/cdip/realtime/';
nameurl = [cdipname 'p1_rt.nc.ascii?'];

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%% FIRST: get tinfo and limit time series if need be
[~, ~, starti, endi, Nt] = get_tinfo(baseurl, cdipname, tlims,mode);
% check if starti, endi, tres, are even/logical
tt = [starti:tres:endi];
% endi = tt(end);
tstring = [ '[' num2str(starti) ':' num2str(tres) ':' num2str(endi) ']'];
tstring = [ '%5B' num2str(starti) ':' num2str(tres) ':' num2str(endi) '%5D']; 

Nt = Nt/tres+1;
if endi~=starti & endi>starti
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% SECOND: get frequency

    % Load frequency
    paramurl = ['waveFrequency' fstr];
    url = [baseurl nameurl paramurl]; %full url name
    f_cdip = dload(url);

    % Load bandwidthh
    paramurl = ['waveBandwidth' fstr];
    url = [baseurl nameurl paramurl]; %full url name
    bw_cdip = dload(url);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% AND THEN: create empty variable, url name, load data, reshape

    %%%% waveEnergyDensity %%%
    if dim==2 
        paramurl = ['wave' varname tstring fstr];
        var_cdip = NaN(fN,Nt);
        url = [baseurl nameurl paramurl]; %full url name
        data = urlread(url); data = strsplit(data,'\n');
        dstart = find(strcmp(data,'---------------------------------------------'))+1;

        for j=1:Nt
            waveVar = strsplit(data{dstart+j}, ', '); 
            if vo; disp(waveVar{1}); end
            waveVar = cellfun(@str2double,waveVar);
            try
                var_cdip(:,j) = waveVar(2:end);
            catch
                1+1
            end
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    spec_cdip_2 = var_cdip;

    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% THIRD: provide time variable
    tstart = find(strcmp(data,['wave' varname '.waveTime[' num2str(Nt) ']']));
    waveVar = strsplit(data{tstart+1}, ', ');
    waveVar = cellfun(@str2double,waveVar);
    toff = datenum(1970,1,1,0,0,0);
    t_cdip_2 = waveVar./(24*60*60) + toff;
else
    spec_cdip_2 = [];
    t_cdip_2 = [];
    disp('---> Not using CDIP realtime');
end

%% COMBINE REALTIME AND HISTORIC
% -------------------------------------------------------------------------
try
    spec_cdip = [spec_cdip_1 spec_cdip_2];
    t_cdip = [t_cdip_1 t_cdip_2];

catch
    spec_cdip = [spec_cdip_2];
    t_cdip = [t_cdip_2];

end


end

function [tstart tend starti endi Nt] = get_tinfo(baseurl, cdipname, tlims,mode)
    if strcmp(mode, 'hist')
        infourl = [baseurl cdipname 'p1/' cdipname 'p1_historic.nc.html'];
%         ilines = [122, 388,389];
    elseif strcmp(mode, 'rt')
        infourl = [baseurl cdipname 'p1_rt.nc.html'];
%         ilines = [120, 398,399];
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
    idx = find(contains(info, 'Int32 waveTime[waveTime = ')); idx = idx(1);
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
    yr = str2num(tstr(1:4)); mo = str2num(tstr(6:7)); day = str2num(tstr(9:10)); hr = str2num(tstr(12:13));
    tend = datenum(yr,mo,day, hr, 0,0);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%     disp(['     CDIP runs from ' datestr(tstart) ' to '  datestr(tend)])
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %%% FIND VAR = starti, endi
    if ~isempty(tlims)
        endi = Nt;
        sdi = 100; % change this if you want to make the time cuttoff more accurate
        starti = Nt-floor(Nt/sdi)*sdi;
        % downloading coarse resolution time series data to get approximate index location for tlims
        % downloading coarse resolution time series data to get approximate index location for tlims
        [simpletime, idx] = get_simpletime(baseurl, cdipname, mode, starti, endi, sdi);
        
        [~, idxstart] = min(abs(simpletime-tlims(1)));
        [~, idxend] = min(abs(simpletime-tlims(2)));
        starti = idx(idxstart); endi = idx(idxend);
        
        if starti~=endi  
            di = 50;
            % fine tune starti
            mi = starti;
            sti = mi - di;  if sti<1; sti = 1; end
            eni = mi + di;  if eni>Nt; eni = Nt; end
            sdi = 1;
            [simpletime, idx] = get_simpletime(baseurl, cdipname, mode, sti, eni, sdi);           
            [~, idxstart] = min(abs(simpletime-tlims(1)));
            starti = idx(idxstart); 
            
            % fine tune endi
            mi = endi;
            sti = mi - di;  if sti<1; sti = 1; end
            eni = mi + di;  if eni>Nt; eni = Nt; end
            sdi = 1;
            [simpletime, idx] = get_simpletime(baseurl, cdipname, mode, sti, eni, sdi);           
            [~, idxend] = min(abs(simpletime-tlims(2)));
            endi = idx(idxend); 
            
        end
        Nt = endi-starti;
        
    else % or small subset (data way to big to pull entire thing)
        starti = Nt-1000;
        endi = Nt;
    end
end


function [simpletime, idx] = get_simpletime(baseurl, cdipname, mode, starti, endi, sdi)
    tvarstr = 'waveTime';
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
end

function [var] = dload(url) % download thredds data from url
    data = urlread(url); data = strsplit(data,'\n');
    dstart = find(strcmp(data,'---------------------------------------------'))+1;

    var = strsplit(data{dstart+1}, ', ');
    var = cellfun(@str2double,var);
end