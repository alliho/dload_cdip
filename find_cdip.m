function [cdipids lons lats] = find_cdip(lonlims, latlims, varargin)
    

%% parse variable inputs
p = inputParser;
addParameter(p, 'archive', 1);
addParameter(p, 'realtime', 1);
parse(p, varargin{:});
rltm = p.Results.realtime;
archv = p.Results.archive;

%% main function

if archv
    rawdat = urlread('https://cdip.ucsd.edu/offline/wavecdf/wnc_listing.php?ARCHIVE');
    dat = strsplit(rawdat, '\n');
    [ii] = find(contains(dat,'p1')); 
    datids = dat(ii);
    datids = cellfun(@(s) [extractBetween(s,'ARCHIVE/', 'p1">')],datids);
    cdipids_archv = cellfun(@str2num, datids);
else
    cdipids_archv = [];
end
if rltm
    rawdat = urlread('https://cdip.ucsd.edu/offline/wavecdf/wnc_listing.php?REALTIME');
    dat = strsplit(rawdat, '\n');
    [ii] = find(contains(dat,'p1_rt')); 
    datids = dat(ii);
    datids = cellfun(@(s) [extractBetween(s,'REALTIME/', 'p1_rt')],datids);
    cdipids_rltm = cellfun(@str2num, datids);
else
    cdipids_rltm = [];
end

cdipids = unique(union(cdipids_archv, cdipids_rltm));

remids = [230 253 268];
cdipids = setdiff(cdipids, remids);
cdipids = sort(unique(cdipids));


%%% SUBSET BASED ON REGION
[lats lons] = arrayfun(@(id) dload_cdiplatlon(id), cdipids);
ii = find(isin(lons, lonlims) & isin(lats, latlims));
cdipids = cdipids(ii); lons = lons(ii); lats = lats(ii);
   

end