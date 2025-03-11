% downloading realtime CDIP data EX
% this function walks through three different examples of downloading CDIP
% data using the dload_cdip package. 
%    1) dload_cdipbuoy  |   download entire buoy with all variables, 
%                           omitting specified variables
%    2) dload_cdipvar   |   download individual variables
%    3) dload_cdipxyz   |   download raw displacement data
%    4) dload_ww3var    |   download ww3 model output at buoy 
%                           (CAUTION: CODE NEEDS TO BE UPDATED)

%% first, add path
addpath(genpath('./dload_cdip_v4/'))

%% EX 1): dload_cdipbuoy
cdipid = 067;
daterange = now + [-20 0];
tres = 1;

% downloads all parameters, omiting fourier coefficients ------------------
clear cdip;
cdip = dload_cdipbuoy(cdipid, daterange,tres, 'omit', {'a1', 'b1', 'a2', 'b2', 'check', 'md'});

% plot --------------------------------------------------------------------

figure(1); clf; 
subplot(211);
plot(cdip.time, cdip.hs, 'k.')
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
title([datestr(cdip.time(1)) ' to ' datestr(cdip.time(end)) ' [UTC]'])
grid on; box on;

subplot(212);
pcolor(cdip.time, cdip.f, log10(cdip.sf)); shading flat;
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
grid on; box on; set(gca, 'Layer', 'top')

%% EX 2): dload_cdipvar (with currents)
cdipid = 132;
daterange = now + [-20 0];
tres = 1;

% spectral data etc. ------------------------------------------------------
clear cdip;

varnm = 'Time';
[cdip.time] = dload_cdipvar(cdipid,varnm,daterange, tres);
varnm = 'Hs';
[cdip.hs] = dload_cdipvar(cdipid,varnm,daterange, tres);

varnm = 'Tp';
[cdip.tp] = dload_cdipvar(cdipid,varnm,daterange, tres);

varnm = 'gpsLongitude';
[cdip.lon] = dload_cdipvar(cdipid,varnm,daterange, tres);

varnm = 'gpsLatitude';
[cdip.lat] = dload_cdipvar(cdipid,varnm,daterange, tres);

varnm = 'acmTime';
[cdip.currents.time] = dload_cdipvar(cdipid,varnm,daterange, tres);

varnm = 'acmDirection';
[cdip.currents.direction] = dload_cdipvar(cdipid,varnm,daterange, tres);

varnm = 'acmSpeed';
[cdip.currents.speed] = dload_cdipvar(cdipid,varnm,daterange, tres);

% plot --------------------------------------------------------------------

figure(2); clf; 
subplot(311);
plot(cdip.time, cdip.hs, 'k.')
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
title([datestr(cdip.time(1)) ' to ' datestr(cdip.time(end)) ' [UTC]'])
grid on; box on;

subplot(312);
plot(cdip.time, cdip.tp, 'k.')
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
grid on; box on;

subplot(313);
plot(cdip.currents.time, cdip.currents.speed, 'k.')
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
grid on; box on;

%% EX 3): dload_cdipxyz
cdipid = 67;
tres = 1;
daterange = [datenum(2021,4,8) datenum(2021, 4,10)]; 

% spectral data etc. ------------------------------------------------------
clear cdip;

varnm = 'Time';
[cdip.time] = dload_cdipvar(cdipid,varnm,daterange, tres);
varnm = 'Hs';
[cdip.hs] = dload_cdipvar(cdipid,varnm,daterange, tres);

% displacement data -------------------------------------------------------
clear cdipxyz;
varname = 'ZDisplacement';
[cdipxyz.(varname), cdipxyz.time] = dload_cdipxyz(cdipid,varname,daterange, tres);
varname = 'XDisplacement';
[cdipxyz.(varname), cdipxyz.time] = dload_cdipxyz(cdipid,varname,daterange, tres);
varname = 'YDisplacement';
[cdipxyz.(varname), cdipxyz.time] = dload_cdipxyz(cdipid,varname,daterange, tres);

% plot --------------------------------------------------------------------

figure(3); clf; 

subplot(3,1,1);hold on;
plot(cdip.time, cdip.hs, 'k.')
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
grid on; box on;
title([datestr(cdip.time(1)) ' to ' datestr(cdip.time(end)) ' [UTC]'])

subplot(3,1,2:3);hold on;
plot(cdipxyz.time, cdipxyz.ZDisplacement, 'k-')
plot(cdipxyz.time, cdipxyz.XDisplacement, 'r-')
plot(cdipxyz.time, cdipxyz.YDisplacement, 'b-')
datetick(gca, 'x')
grid on; box on;
xlim(cdip.time([1 end]))

%% EX 4): dload_ww3var 
cdipid = 067;
daterange = now + [-20 0];
tres = 1;

% CDIP data ---------------------------------------------------------------
clear cdip;

varnm = 'Time';
[cdip.time] = dload_cdipvar(cdipid,varnm,daterange, tres);
varnm = 'Hs';
[cdip.hs] = dload_cdipvar(cdipid,varnm,daterange, tres);

% WW3 data ---------------------------------------------------------------
clear ww3;
ndbcid = get_ndbcid(cdipid);

varnm = 'Time';
[ww3.time] = dload_ww3var(ndbcid,varnm,daterange, tres);
varnm = 'Hs';
[ww3.hs] = dload_ww3var(ndbcid,varnm,daterange, tres);

% plot --------------------------------------------------------------------
figure(4); clf; 

subplot(1,1,1);hold on;
plot(cdip.time, cdip.hs, 'k.', 'DisplayName','CDIP')
plot(ww3.time, ww3.hs, 'r.', 'DisplayName','WW3')
xlim(cdip.time([1 end]))
datetick(gca, 'x', 'keepticks', 'keeplimits')
grid on; box on;
legend
title([datestr(cdip.time(1)) ' to ' datestr(cdip.time(end)) ' [UTC]'])



