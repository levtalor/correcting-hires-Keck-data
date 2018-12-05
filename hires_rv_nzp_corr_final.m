function [stars,Nrv,rvc_mean_vec,rvc_mean_vec_err,Nnoaverage,Nout_star,Nout_night,ind_same_vec,rv_std_vec,rv_med_err,nights_mean,nights_mean_err,nights_mean_std,Nrv_night,rv_std_vec_corr,rv_med_err_corr,nzp,enzp] = ...
    hires_rv_nzp_corr_final(fname,jd_min,jd_max,jd_jump,Nrv_min,Nrv_min_night,sigma_outlier_star,sigma_outlier_night,init_rv_var,final_rv_var,init_rv_systematics,jump_corr,smooth_nzp,smooth_span,...
                         drift_corr,showflag,saveflag,savedata,dirname,save_fmt,unselfbias)
% [stars,Nrv,rvc_mean_vec,rvc_mean_vec_err,Nnoaverage,Nout_star,Nout_night,ind_same_vec,rv_std_vec,rv_med_err,nights_mean,nights_mean_err,nights_mean_std,Nrv_night,rv_std_vec_corr,rv_med_err_corr,nzp,enzp] = ...
%    hires_rv_nzp_corr_02(fname,jd_min,jd_max,jd_jump,Nrv_min,Nrv_min_night,sigma_outlier_star,sigma_outlier_night,init_rv_var,final_rv_var,init_rv_systematics,jump_corr,smooth_nzp,smooth_span,...
%                         drift_corr,showflag,saveflag,savedata,dirname,save_fmt,unselfbias)
% Corrects for small systematic effects in HIRES RVs published by Butler+17.
%
% Input:
% fname = full path to the file that contains full pathes to the hires RV
%         files (default = 'C:\lev\HIRES\vel_files.txt').
% jd_min - first night of real survey (default = 2450275 = Jul 10, 1996)
% jd_max - last night of public data (default = 2457004 = March 03, 2014)
% jd_jump = 2453225 - night of HIRES upgrade (Aug 07, 2004)
% Nrv_min - minimum RVs/star to be included in the analysis (default = 5).
% Nrv_min_night - minimum RVs/night to correct for its average (default = 3).
% sigma_outlier_star - Nsigma for outlier rejection per star (default = 10).
% sigma_outlier_night - Nsigma for outlier rejection per night (default = 5).
% init_rv_var - initial minimum std(RV) of an RV-loud star (m/s, default = 10).
% final_rv_var - final minimum std(RV) of an RV-loud star (m/s, default = 10).
% init_rv_systematics - initial guess of the systematic RV uncertainty
%                       (added to RV uncertainties, default = 1.0 m/s).
% jump_corr - flag to correct Aug 2004 upgrade as jump in nzps (default = 0).
% smooth_nzp - flag to smooth NZPs and use the smothed version as NZPs (default = 1).
% smooth_span - smoothing span (default = 50 days).
% drift_corr - flag to correct glogal nightly drift (default = 1).
% showflag - flag to make plots (default = 1).
% saveflag - flag to save plots (default = 1).
% savedata - flag to save data (default = 1).
% dirname - directory to store the results (default = 'C:\lev\HIRES\results\').
% save_fmt - format of the plots (default = '.fig').
% unselfbias - flag to avoid self biasing (default = 0).
%
% Output:
% stars - cell vector with the names of the stars
% Nrv - number of measurements per star
% rvc_mean_vec - vector of mean RVC per star
% rvc_mean_vec_err - vector of the error of the mean RVC per star
% Nnoaverage - number of measurements with no nightly-average correction
% Nout_star - Number of outliers per-star
% Nout_night - Number of outliers per-night
% ind_same_vec - number of nights in which the star was observed more than once
% rv_std_vec - std(RV) per star (before correction)
% rv_med_err - median RV error per star (before correction)
% nights_mean - nightly means
% nights_mean_err - nightly-mean error
% nights_mean_std - nightly-mean scatter divided by sqrt(Nrv_night)
% Nrv_night - number of good RVs per night
% rv_std_vec_corr - std(RV) per star (after correction)
% rv_med_err_corr - median RV error per star (after correction)
% nzp - NZPs before smoothing (including NaNs for no-mean nights)
% enzp - NZP errors before smoothing (including NaNs for no-mean nights)
%
% Created: 20160609 LT
% Last modified: 20181114 LT
%

%%
tic

%% input parameters:
if nargin<1
    fname = 'C:\lev\HIRES\vel_files.txt'; 
    sigma_outlier_star = 10; % Nsigma for outlier rejection per star (default = 10).
    sigma_outlier_night = 5; % Nsigma for outlier rejection per night (default = 5).
    Nrv_min = 5; % minimum RVs/star for it to be included in the analysis (default = 5).
    Nrv_min_night = 3; % minimum RVs/night to correct for its average (default = 3).
    jd_min = 2450275; % first night of the public survey data (July 10, 1996)
    jd_max = 2457004; % last night of the public survey data (March 03, 2014)
    jd_jump = 2453225; % night of HIRES upgrade (Aug 07, 2004)
    init_rv_var = 10; % initial rv std of a variable star (m/s), (default = 10).
    final_rv_var = 10; % final rv std of a variable star (m/s), (default = 10).
    init_rv_systematics = 1.0; % initial guess of the systematic rv uncertainty (added to RV uncertainties).
    jump_corr = 0; % flag to correct Aug 2004 upgrade as jump in nzps.
    smooth_nzp = 1; % flag to smooth NZPs and use the smoothed version as NZPs.
    smooth_span = 50; % smoothing span (days)
    drift_corr = 1; % flag to correct for the average nightly drift
    showflag = 0; % flag to make plots
    saveflag = 0; % flag to save plots
    savedata = 1; % flag to save data
    dirname = 'C:\lev\HIRES\results\'; % directory to store the results
    save_fmt = '.fig';
    unselfbias = 1; % flag to avoid self biasing (NOTE: runs for ~2.4 hr and makes no plots + disable remove_var option).
    remove_var = 0; % flag to remove specific stars to check self biasing
    variables = {'HD7924'};
end
if remove_var
    dirname = [dirname 'exclude' variables{1} '\'];
end

%% ===========================================================================
% STARTING POINT OF CALCULATING THE SPARSE MATRIX OF RV-NIGHT
%===========================================================================

% load the list of .vel files and initialize some output vectors:
targets = loadtxt(fname);
target_num = length(targets);
stars = loadtxt('stars.txt');
rv_std_vec = nan(length(targets),1);
Nrv = rv_std_vec;
rvc_mean_vec = Nrv;
rvc_mean_vec_err = Nrv;
rv_med_err = Nrv;
ind_same_vec = Nrv;
Nnoaverage = zeros(length(targets),1);
Nout_star = Nrv;

% Initialize the matrices:
a = clock;
day_str = [num2str(a(1),'%04d') num2str(a(2),'%02d') num2str(a(3),'%02d')];
% jd_now = utc2jd(a(1),a(2),a(3),a(4),a(5),a(6));
m = jd_max-jd_min+1;
rv_mat = nan(target_num,m);
err_mat = rv_mat;
jd_mat = rv_mat;

% open a .orb file:
if savedata
    orbname = [dirname 'orb_files\hires_preNZP.orb'];
    fid = fopen(orbname,'w');
end

%% collect the times, RVs, and errors into the matrices, target by target:
for t=1:target_num
    
    % load the .vel file:
    filename = targets{t};
    M = load(filename);
    
    % derive the starname:
    starname = stars{t};
    
    % Empty files could be caused by not-drift-corrected measurement:
    if isempty(M)
        warning([starname '.vel is an empty file']);
        continue
    end
    
    % read matrix into vectors:
    M = sortrows(M,1);
    bjd_tmp = M(:,1)-0.5;
    M(bjd_tmp<jd_min,:)=[]; % Remove measurements prior to jd_min:
    M(bjd_tmp>jd_max,:)=[]; % Remove measurements after jd_max:
    bjd = M(:,1)-0.5; % Subtract 0.5 day from time stamps to have all 
                      % the RVs of each night with the same floor(bjd).
    rvc = M(:,2);
    rvc_err = M(:,3);
    rv_med_err(t) = nanmedian(rvc_err);

    % Find outliers (do not rely on RV errors - they do not account for systematics yet):
    rv_median = nanmedian(rvc);
    bias = 1-1/4/length(rvc);
    rv_std = 1.48*mad(rvc,1)/bias;
    bad_rv = find(abs(rvc-rv_median)>sigma_outlier_star*rv_std & abs(rvc-rv_median)>init_rv_var);
    rv_std_vec(t) = rv_std;
    Nrv(t) = length(rvc);
    Nout_star(t) = length(bad_rv);
            
    % Write an entry to .orb file (original RVs including the outliers):
    if savedata
        fprintf(fid,'STAR: %s\n',starname);
        for i=1:length(bjd)
            fprintf(fid,'A %f %f %f\n',bjd(i)+0.5,0.001*rvc(i),0.001*rvc_err(i));
        end
        fprintf(fid,'END\n');
    end
    
    % Remove stellar outliers (only for the purpose of calculating the NZPs):
    bjd(bad_rv) = [];
    rvc(bad_rv) = [];
    rvc_err(bad_rv) = [];
    
    % For nights with multiple exposures, take the weighted mean RVC in that
    % night, with its scatter as an error
    nights = floor(bjd);
    columns = nights-jd_min+1;
    ind_same_tmp = 0;
    for i=1:m
        ind_same = find(abs(columns-i)==0);
        if ~isempty(ind_same) && length(ind_same)>1
            [~,ind_min_tmp] = min(rvc_err(ind_same));
            % calculate a weihted mean values:
            bjd_ind_same = nanwmean(bjd(ind_same),rvc_err(ind_same).^-2);
            rvc_ind_same = nanwmean(rvc(ind_same),rvc_err(ind_same).^-2);
            rvc_err_ind_same = max(rvc_err(ind_same(ind_min_tmp)),wstd(rvc(ind_same),rvc_err(ind_same).^-2));
            % save the corrected values instead of the minimum err value:
            bjd(ind_same(ind_min_tmp)) = bjd_ind_same;
            rvc(ind_same(ind_min_tmp)) = rvc_ind_same;
            rvc_err(ind_same(ind_min_tmp)) = rvc_err_ind_same;
            % erase the index from the values to be erased:
            ind_same(ind_min_tmp) = [];
            % erase the the rest of the values:
            columns(ind_same) = [];
            nights(ind_same) = [];
            bjd(ind_same) = [];
            rvc(ind_same) = [];
            rvc_err(ind_same) = [];
            % count the number of binned data
            ind_same_tmp = ind_same_tmp+length(ind_same);
        end
    end
    ind_same_vec(t) = ind_same_tmp;
    
    % Re-calculate the weighted mean (zero-point of the star) and error:
    if Nrv(t)>0
        bias = 1-1/4/length(rvc);
        rvc_err = sqrt(rvc_err.^2+init_rv_systematics^2);
        weights = rvc_err.^-2;
        rvc_mean = nanwmean(rvc,weights);
        rvc_mean_vec(t) = rvc_mean;
        rvc_mean_err = 1/sqrt(sum(weights));
        rvc_mean_std = wstd(rvc,weights)/sqrt(Nrv(t))/bias;
        rvc_mean_err = max(rvc_mean_err,rvc_mean_std);
        rvc_mean_vec_err(t) = rvc_mean_err;

        % Write into the matrices:
        jd_mat(t,columns) = bjd;
        rv_mat(t,columns) = rvc-rvc_mean; % stellar zero-point subtracted RVs
        err_mat(t,columns) = sqrt(rvc_err.^2+rvc_mean_err^2); % the error on the stellar zero-point is co-added to the RV error,
                                                              % so that variables and/or faint stars get lower weight in NZP calculation. 
    else
        warning([starname ' has no good RVs.']);
    end
    
end

% save the matrices
save([dirname 'mat_files\jd_mat'],'jd_mat','rv_mat','err_mat');

% close the .orb and all_rvs files
if savedata
    fclose(fid);
end

%% =========================================================================
% STARTING POINT OF CORRECTING THE RADIAL VELOCITIES WITH SELF BIASING
%===========================================================================

if ~unselfbias
% Find variable stars:
    novar_rv = find(rv_std_vec<init_rv_var);
    rv_std_median = nanmedian(rv_std_vec(novar_rv));
    if showflag
        disp('Before correction:');
        fprintf(1,'Median std(RV) of RV-quiet stars: %4.2f\n',rv_std_median);
        fprintf(1,'std(RV) threshold for an RV_loud star: %4.2f\n',init_rv_var);
    end
    var_rv = find(rv_std_vec>init_rv_var);
    
    % Remove also hidden-variable stars:
    for t=1:target_num
        starname = stars{t}; 
        for v = 1:length(variables)
            if strcmp(variables{v},starname) && remove_var
               var_rv = [var_rv;t];
            end
        end
    end
    
    % remove multiple entries:
    var_rv = sort(var_rv);
    var_rv(diff(var_rv)==0)=[];
    
    
%% plot std(RV) histogram
if showflag
     figure(16);
     set(16,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
     hold off
     g = logspace(0,1.05*log10(max(rv_std_vec)),ceil(10*(1.05*log10(max(rv_std_vec)))));
     hb = hist(rv_std_vec(Nrv>=Nrv_min),g);
     stairs(g,hb,'r','LineWidth',2);
     hold on;
     xlim([1 10^(1.05*log10(max(rv_std_vec)))]);
     ylim([0 1.05*max(hb)]);
     set(gca,'Xscale','log');
     set(gca,'xTickLabel',[1 10 100 1000 10000]);
     ylabel('Number of stars');
     xlabel('std(RV) [m/s]');
     title([date ': Histograms of std(RV) before and after correcting for the NZPs']);
     set(gca,'FontSize',16);
     grid on
end

%% Average RVs from the same night using all stars observed in that night:

% initialize some vectors:
nights_mean = nan(1,m);
nights_median = nights_mean;
nights_mean_err = nights_mean;
nights_mean_std = nights_mean;
Nrv_night = nan(1,m);
Nout_night = Nrv_night;

% load the results of stage 1 (the sparse RV-night matrices):
load([dirname 'mat_files\jd_mat.mat'],'jd_mat','rv_mat','err_mat');

% Remove the RVs of targets with less than Nrv_min exposures (only for the purpose of calculating the NZPs):
jd_mat(Nrv<Nrv_min,:) = nan;
rv_mat(Nrv<Nrv_min,:) = nan;
err_mat(Nrv<Nrv_min,:) = nan;

% Remove the RVs of variable targets (only for the purpose of calculating the NZPs):
jd_mat(var_rv,:) = nan;
rv_mat(var_rv,:) = nan;
err_mat(var_rv,:) = nan;  

% start the plot of all RVs of RV-quiet stars and the NZPs:
if showflag
    figure(10);
    set(10,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    hold off
end
for t=1:m
    night_jds = jd_mat(:,t);
    night_rvs = rv_mat(:,t);
    night_err = err_mat(:,t);
    
    % plot nightly RVs
    if showflag
        figure(10);
        hold on
        plot(night_jds-jd_min+0.5,night_rvs,'.c');
    end
    
    % remove outliers per night (but plot them) only if there are enough RVs:
    good_rvs_tmp = find(isfinite(night_rvs));
    Nrv_night_tmp = length(good_rvs_tmp);
    if Nrv_night_tmp>Nrv_min_night
        rv_median_night = nanmedian(night_rvs);
        bias = 1-1/4/Nrv_night_tmp;
        rv_std_night = 1.48*mad(night_rvs,1)/bias;
        bad_rv_night = find(abs(night_rvs-rv_median_night)>sigma_outlier_night*rv_std_night & abs(night_rvs-rv_median_night)>init_rv_var);
        if showflag
            plot(night_jds(bad_rv_night)-jd_min+0.5,night_rvs(bad_rv_night),'or')
        end
        night_jds(bad_rv_night) = nan;
        night_rvs(bad_rv_night) = nan;
        night_err(bad_rv_night) = nan;
        Nout_night(t) = length(bad_rv_night);
    end
    
    %==========================================================
    % Calculate the nightly zero-point RV using only good RVs:
    good_rvs = find(isfinite(night_rvs));
    Nrv_night(t) = length(good_rvs);
    bias = 1-1/4/Nrv_night(t);
    weights = night_err.^-2;
    if sum(isfinite(weights))
        nights_mean(t) = nanwmean(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_err(t) = 1/sqrt(sum(weights(good_rvs)))/bias;
        nights_median(t) = wmedian(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_std(t) = wstd(night_rvs(good_rvs),weights(good_rvs))/sqrt(Nrv_night(t))/bias;
        % Take the max of nights_mean_err and nights_mean_std:
        nights_mean_err(t) = max(nights_mean_err(t),nights_mean_std(t));
    end
    %==========================================================
    
end

%% display some statistics:
    median_nights_mean = nanmedian(nights_mean(Nrv_night>=Nrv_min_night))
    std_nights_mean = 1.48*mad(nights_mean(Nrv_night>=Nrv_min_night),1)

%% find nights where NZP correction could not be applied and mark them:
bad_mean = find(Nrv_night<Nrv_min_night);
nights_vec = 1:m;

if showflag
    figure(10)
    h = errorbar(nights_vec(isfinite(nights_mean)),nights_mean(isfinite(nights_mean)),nights_mean_err(isfinite(nights_mean)),'.k','MarkerSize',15,'LineWidth',1);
% Make cap-lines mult times long
    mult = 0;                               % elongation factor
    b = h.Bar;                              % hidden property of h=errorbar(X,Y,E)
    drawnow                                 % populate b's properties
    vd = b.VertexData;
    X = nights_vec(isfinite(nights_mean));
    N = numel(X);                           % number of error bars
    capLength = vd(1,2*N+2,1) - vd(1,1,1);  % assumes equal length on all
    newLength = capLength * mult;
    leftInds = N*2+1:2:N*6;
    rightInds = N*2+2:2:N*6;
    vd(1,leftInds,1) = [X-newLength, X-newLength];
    vd(1,rightInds,1) = [X+newLength, X+newLength];
    b.VertexData = vd;
% plot the bad-mean NZPs:
    plot(nights_vec(bad_mean),nights_mean(bad_mean),'sr');
end

%% ==========================================================
% Put NaN for no-mean nights for the smoothing
nights_mean(bad_mean) = NaN;
nights_mean_err(bad_mean) = NaN;

% smooth the NZPs:
    nzp = nights_mean(:);
    enzp = nights_mean_err(:);
    weight = 1./enzp.^2;
    if ~jump_corr && smooth_nzp
        [nzp_filt,nzp_filt_err] = wmeanfilt02(nights_vec,nzp,weight,smooth_span);
        nights_mean = nzp_filt';
        nights_mean_err = nzp_filt_err';
        % linearly interpolate between nights of no NZP:
        nights_mean = interp1(nights_vec(~isnan(nzp)),nights_mean(~isnan(nzp)),nights_vec,'linear',0);
        nights_mean_err = interp1(nights_vec(~isnan(nzp)),nights_mean_err(~isnan(nzp)),nights_vec,'linear',1);
        if showflag
            plot(nights_vec,nights_mean,'-g','LineWidth',2.5);
        end
    elseif jump_corr && smooth_nzp
        [nzp_filt_b,nzp_filt_err_b] = wmeanfilt02(nights_vec(1:(jd_jump-jd_min)),nzp(1:(jd_jump-jd_min)),weight(1:(jd_jump-jd_min)),smooth_span);
        [nzp_filt_a,nzp_filt_err_a] = wmeanfilt02(nights_vec(jd_jump-jd_min+1:end),nzp(jd_jump-jd_min+1:end),weight(jd_jump-jd_min+1:end),smooth_span);
        nzp_filt_jump = [nzp_filt_b;nzp_filt_a];
        nzp_filt_err = [nzp_filt_err_b;nzp_filt_err_a];
        nights_mean = nzp_filt_jump';
        nights_mean_err = nzp_filt_err';
        % linearly interpolate between nights of no NZP:
        nights_mean = interp1(nights_vec(~isnan(nzp)),nights_mean(~isnan(nzp)),nights_vec,'linear',0);
        nights_mean_err = interp1(nights_vec(~isnan(nzp)),nights_mean_err(~isnan(nzp)),nights_vec,'linear',1);
        if showflag
            plot(nights_vec,nights_mean,'-g','LineWidth',2);
        end
    else
        % put medin and satter for no-mean nights: 
        nights_mean(bad_mean) = median_nights_mean;
        nights_mean_err(bad_mean) = std_nights_mean;
    end
% ==========================================================
    
%% some design:
if showflag
    figure(10);
    xlim([-100 6850]);
    hold on;
    plot([-100 6850],[10 10],'-k');
    ylim([-10 10]);
    plot([6850 6850],[-10 10],'-k');
    xlabel (['JD - ' num2str(2450275)],'fontsize',14);
    set(10,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    set(gca,'FontSize',12)
    ylabel('RV [m/s]');
    set(gca,'FontSize',12);
    grid on
    if saveflag
        filename = ['hires_night_means_selfbiased_' day_str];
        saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end
end

%% calculate the global nightly drift
ObsCoo = -0.431883; % Longitude of the telescope
jd = nan(target_num*m,1);
rva = jd;
err = jd;
rv_corr = rv_mat;
for j=1:m;
    rv_corr(:,j)=rv_corr(:,j)-nights_mean(j);
end
for i=1:target_num;
    jd(m*(i-1)+1:m*(i-1)+m)=jd_mat(i,:)';
    rva(m*(i-1)+1:m*(i-1)+m)=rv_corr(i,:)';
    err(m*(i-1)+1:m*(i-1)+m)=err_mat(i,:)';
end
inan = isnan(jd);
jd(inan)=[];
rva(inan)=[];
err(inan)=[];
clear inan
modbjd = mod(jd,1)+ObsCoo; % jd is bjd-0.5; 
% =====================================================================
if drift_corr
   [a1,b1,da1,db1,~,~,~,~,~,Ftest1] = linfit(modbjd,rva,err,0);
   if showflag
    figure(11)
        plot(modbjd,rva,'.')
        set(11,'units','normalized','outerposition',[0.05 0.05 0.6 0.9]);
        hold on
        grid on
        plot([min(modbjd)-0.01 max(modbjd)+0.01],b1+[min(modbjd)-0.01 max(modbjd)+0.01]*a1,'-k','LineWidth',2.5)
        ylabel('NZP-corrected RV-quiet star RVs [m/s]');
        ylim([-10 10])
        xlim([min(modbjd)-0.01 max(modbjd)+0.01]);
        set(gca,'xTick',[-3/12 -2/12 -1/12 0 1/12 2/12 3/12]);
        set(gca,'xTickLabel',[-6 -4 -2 0 2 4 6]);
        xlabel('t_{mid} [hr]');
        legend(['F_{test} = ' num2str(Ftest1) ' (' num2str(length(rva)) ' points)'],['line fit: a = ' num2str(a1) '+/-' num2str(da1)]);
        set(gca,'FontSize',16)
        if saveflag
           filename = ['hires_nightly_drift_selfbiased_' day_str];
           saveas(gcf,[dirname 'figures\' filename save_fmt]);
        end
   end
end
% =====================================================================
%% apply the nightly correction, add its uncertainties, and write a new .orb and avc.vel files:
if savedata
    orbname = [dirname 'orb_files\hires_avc_selfbiased_' day_str '.orb'];
    fid = fopen(orbname,'w');
    allname = [dirname 'txt_files\hires_all_rvs_selfbiased_' day_str '.tsv'];
    fid_all = fopen(allname,'w');
end
rv_std_vec_corr = nan(length(targets),1);
rv_med_err_corr = nan(length(targets),1);

for t=1:target_num
   
    % Re-load all the RVs from the dat file for correction:
    filename = targets{t};
    M = load(filename);
    
    % derive the starname:
    starname = stars{t};
    
    % skip empty files:
    if isempty(M)
        warning([starname '.vel is an empty file']);
        continue
    end
    
    % read matrix into vectors:
    M = sortrows(M,1);
    bjd_tmp = M(:,1)-0.5;
    M(bjd_tmp<jd_min,:)=[]; % Remove measurements prior to jd_min:
    bjd = M(:,1)-0.5;
    modbjd_t = mod(bjd,1)+ObsCoo;
    rvc = M(:,2);
    rvc_err = M(:,3);
    col4 = M(:,4);
    col5 = M(:,5);
    col6 = M(:,6);
    col7 = M(:,7);
    corr_t = nan(size(bjd));
    corr_err = nan(size(bjd));
    
    %==========================================================
    % subtract the nightly mean and co-add its error:
    nights = floor(bjd);
    columns = nights-jd_min+1;
    for i=1:length(bjd)
        if drift_corr
            corr_t(i) = nights_mean(columns(i)) + (b1+modbjd_t(i)*a1);
        else
            corr_t(i) = nights_mean(columns(i));
        end
        rvc(i) = rvc(i) - corr_t(i);
        if drift_corr
            corr_err(i) = sqrt(nights_mean_err(columns(i))^2 + db1^2 + (modbjd_t(i)*da1)^2);
        else
            corr_err(i) = nights_mean_err(columns(i));
        end
        rvc_err(i) = sqrt(rvc_err(i)^2 + corr_err(i)^2);
    end
    %==========================================================
    
    % calculate the std of corrected RVs:
    bias = 1-1/4/length(rvc);
    rv_std = 1.48*mad(rvc,1)/bias;
    rv_std_vec_corr(t) = rv_std;
    rv_med_err_corr(t) = nanmedian(rvc_err);
    Nrv(t) = length(rvc);

   if savedata
    % Write an entry to .orb file (corrected RVs, NO outliers removed):
    fprintf(fid,'STAR: %s\n',starname);
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid,'A %f %f %f\n',bjd(i)+0.5,0.001*rvc(i),0.001*rvc_err(i));
        end
    end
    fprintf(fid,'END\n');
    
    % Write an entry to all_rvs file (corrected RVs, NO outliers removed):
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid_all,'%14s\t%13.5f\t%9.3f\t%7.3f\t%7.4f\t%6.4f\t%6d\t%4d\t%5.3f\t%5.3f\n',starname,bjd(i)+0.5,rvc(i),rvc_err(i),col4(i),col5(i),col6(i),col7(i),corr_t(i),corr_err(i));
        end
    end
    
    % Write a .avc.vel file (corrected RVs, NO outliers removed):
    fid_avc = fopen([dirname 'corrected_rvs/' starname '_selfbiased.avc.vel'],'w');
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid_avc,'%f %f %f\n',bjd(i)+0.5,rvc(i),rvc_err(i));
        else
            Nnoaverage(t) = Nnoaverage(t)+1;
        end
    end
    fclose(fid_avc);
   end

end

if savedata
    fclose(fid);
    fclose(fid_all);
end

%% Find again RV-loud stars and make some display:
rv_std_median_after = nanmedian(rv_std_vec_corr(novar_rv));
if showflag
    disp('After correction:');
    fprintf(1,'Median std(RV) of RV-quiet stars: %4.2f\n',rv_std_median_after);
    fprintf(1,'std(RV) threshold for an RV_loud star: %4.1f\n',final_rv_var);
    fprintf(1,'The actually used RV-std threshold for an RV_loud star: %3.1f\n',final_rv_var);
end    

%% Plot the histogram of std(RV) after correction
if showflag
     figure(16);
     hold on
     g = logspace(0,1.05*log10(max(rv_std_vec_corr)),ceil(10*(1.05*log10(max(rv_std_vec_corr)))));
     ha = hist(rv_std_vec_corr(Nrv>=Nrv_min),g);
     stairs(g,ha,'b','LineWidth',2);
     xlim([1 10^(1.05*log10(max(rv_std_vec_corr)))]);
     ylim([0 1.05*max([ha(:);hb(:)])]);
     set(gca,'Xscale','log');
     set(gca,'xTickLabel',[1 10 100 1000 10000]);
      ylabel('Number of stars');
     legend(['Median(std(RV)) of RV-quiet stars before correction:' num2str(rv_std_median,3) ' m/s'],...
            ['Median(std(RV)) of RV-quiet stars after correction:' num2str(rv_std_median_after,3) ' m/s']);
    if saveflag
        filename = ['hires_stdRV_hist_selfbiased_' day_str];
        saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end   
end

%% plot std(RV) before and after correction for RV-quiet stars:
if showflag
    figure(17)
    plot(rv_std_vec,rv_std_vec_corr,'.b','MarkerSize',15);
    hold on
    ylim([0 final_rv_var])
    xlim([0 final_rv_var])
    plot(0:1:final_rv_var,0:1:final_rv_var,'-k','LineWidth',2);
    [a,b] = linfit(rv_std_vec(novar_rv),rv_std_vec_corr(novar_rv),rv_std_vec_corr(novar_rv)./sqrt(2.*Nrv(novar_rv)));
    plot(rv_std_vec(novar_rv),b+a*rv_std_vec(novar_rv),'-g','LineWidth',2);
    xlabel('std(RV) before correction [m/s]');
    ylabel('std(RV) after correction [m/s]');
    grid on
    set(17,'units','normalized','outerposition',[0.25 0.05 0.6 0.9]);
    set(gca,'FontSize',16);
    title([date ': NZP correction of hires RVs']);
    if saveflag
         filename = ['hires_RVstd_compare_selfbiased_' day_str];
         saveas(gcf,[dirname 'figures\' filename save_fmt]);
    end 
end

%% ===========================================================================
% STARTING POINT OF CORRECTING THE RADIAL VELOCITIES WITHOUT SELF BIASING
%===========================================================================
else % = if unselfbias

if savedata
    orbname = [dirname 'orb_files\hires_avc_' day_str '.orb'];
    fid = fopen(orbname,'w');
    allname = [dirname 'txt_files\hires_all_rvs_' day_str '.tsv'];
    fid_all = fopen(allname,'w');
end

rv_std_vec_corr = nan(length(targets),1);
rv_med_err_corr = nan(length(targets),1);

for t=1:target_num
    % derive the starname:
    starname = stars{t};
    fprintf(1,'Running star number %d/%d (%s)\n',t,target_num,starname);
    
    % Find variable stars:
    var_rv = find(rv_std_vec>init_rv_var);
    
    % Remove also the current star:
    var_rv = [var_rv;t];
    
    % remove multiple entries:
    var_rv = sort(var_rv);
    var_rv(diff(var_rv)==0)=[];
    
%% Average RVs from the same night using all stars observed in that night:

% initialize some vectors:
nights_mean = nan(1,m);
nights_median = nights_mean;
nights_mean_err = nights_mean;
nights_mean_std = nights_mean;
Nrv_night = nan(1,m);
Nout_night = Nrv_night;

% load the results of stage 1 (the sparse RV-night matrices):
load([dirname 'mat_files\jd_mat.mat'],'jd_mat','rv_mat','err_mat');

% Remove the RVs of targets with less than Nrv_min exposures:
jd_mat(Nrv<Nrv_min,:) = nan;
rv_mat(Nrv<Nrv_min,:) = nan;
err_mat(Nrv<Nrv_min,:) = nan;

% Remove the RVs of variable targets (including the current star):
jd_mat(var_rv,:) = nan;
rv_mat(var_rv,:) = nan;
err_mat(var_rv,:) = nan;  

%% =========================================================================
% STARTING POINT OF CALCULATING THE CORRECTION (NESTED LOOP)
%===========================================================================
for y=1:m
    night_jds = jd_mat(:,y);
    night_rvs = rv_mat(:,y);
    night_err = err_mat(:,y);
    
    % remove outliers per night (but plot them) only if there are enough RVs:
    good_rvs_tmp = find(isfinite(night_rvs));
    Nrv_night_tmp = length(good_rvs_tmp);
    if Nrv_night_tmp>Nrv_min_night
        rv_median_night = nanmedian(night_rvs);
        bias = 1-1/4/Nrv_night_tmp;
        rv_std_night = 1.48*mad(night_rvs,1)/bias;
        bad_rv_night = find(abs(night_rvs-rv_median_night)>sigma_outlier_night*rv_std_night & abs(night_rvs-rv_median_night)>init_rv_var);
        night_jds(bad_rv_night) = nan;
        night_rvs(bad_rv_night) = nan;
        night_err(bad_rv_night) = nan;
        Nout_night(y) = length(bad_rv_night);
    end
    
    %==========================================================
    % Calculate the nightly zero-point RV using only good RVs:
    good_rvs = find(isfinite(night_rvs));
    Nrv_night(y) = length(good_rvs);
    bias = 1-1/4/Nrv_night(y);
    weights = night_err.^-2;
    if sum(isfinite(weights))
        nights_mean(y) = nanwmean(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_err(y) = 1/sqrt(sum(weights(good_rvs)))/bias;
        nights_median(y) = wmedian(night_rvs(good_rvs),weights(good_rvs));
        nights_mean_std(y) = wstd(night_rvs(good_rvs),weights(good_rvs))/sqrt(Nrv_night(y))/bias;
        nights_mean_err(y) = max(nights_mean_err(y),nights_mean_std(y));
    end
    %==========================================================
    
end
%% ===========================================================================
% END POINT OF CALCULATING THE CORRECTION (NESTED LOOP)
%===========================================================================

%% display some statistics:
    median_nights_mean = nanmedian(nights_mean)
    std_nights_mean = 1.48*mad(nights_mean,1)
    
%% find nights where NZP correction could not be applied and mark them:
bad_mean = find(Nrv_night<Nrv_min_night);
nights_vec = 1:m;

%% ==========================================================
% Put NaN for no-mean nights for the smoothing
nights_mean(bad_mean) = NaN;
nights_mean_err(bad_mean) = NaN;

% smooth the NZPs:
    nzp = nights_mean(:);
    enzp = nights_mean_err(:);
    weight = 1./enzp.^2;
    if ~jump_corr && smooth_nzp
        [nzp_filt,nzp_filt_err] = wmeanfilt02(nights_vec,nzp,weight,smooth_span);
        nights_mean = nzp_filt';
        nights_mean_err = nzp_filt_err';
        % linearly interpolate between nights of no NZP:
        nights_mean = interp1(nights_vec(~isnan(nzp)),nights_mean(~isnan(nzp)),nights_vec,'linear',0);
        nights_mean_err = interp1(nights_vec(~isnan(nzp)),nights_mean_err(~isnan(nzp)),nights_vec,'linear',1);
    elseif jump_corr && smooth_nzp
        [nzp_filt_b,nzp_filt_err_b] = wmeanfilt02(nights_vec(1:(jd_jump-jd_min)),nzp(1:(jd_jump-jd_min)),weight(1:(jd_jump-jd_min)),smooth_span);
        [nzp_filt_a,nzp_filt_err_a] = wmeanfilt02(nights_vec(jd_jump-jd_min+1:end),nzp(jd_jump-jd_min+1:end),weight(jd_jump-jd_min+1:end),smooth_span);
        nzp_filt_jump = [nzp_filt_b;nzp_filt_a];
        nzp_filt_err = [nzp_filt_err_b;nzp_filt_err_a];
        nights_mean = nzp_filt_jump';
        nights_mean_err = nzp_filt_err';
        % linearly interpolate between nights of no NZP:
        nights_mean = interp1(nights_vec(~isnan(nzp)),nights_mean(~isnan(nzp)),nights_vec,'linear',0);
        nights_mean_err = interp1(nights_vec(~isnan(nzp)),nights_mean_err(~isnan(nzp)),nights_vec,'linear',1);
    else
        % put medin and satter for no-mean nights: 
        nights_mean(bad_mean) = median_nights_mean;
        nights_mean_err(bad_mean) = std_nights_mean;
    end
% ==========================================================

%% calculate the global nightly drift
ObsCoo = -0.431883; % Longitude of the telescope
jd = nan(target_num*m,1);
rva = jd;
err = jd;
rv_corr = rv_mat;
for j=1:m;
    rv_corr(:,j)=rv_corr(:,j)-nights_mean(j);
end
for i=1:target_num;
    jd(m*(i-1)+1:m*(i-1)+m)=jd_mat(i,:)';
    rva(m*(i-1)+1:m*(i-1)+m)=rv_corr(i,:)';
    err(m*(i-1)+1:m*(i-1)+m)=err_mat(i,:)';
end
inan = isnan(jd);
jd(inan)=[];
rva(inan)=[];
err(inan)=[];
clear inan
modbjd = mod(jd,1)+ObsCoo; % jd is bjd-0.5; 
% =====================================================================
if drift_corr
   [a1,b1,da1,db1] = linfit(modbjd,rva,err,0);
end
% =====================================================================
%% apply the nightly correction, add its uncertainties, and write a new .orb and avc.vel files:
   
    % Re-load all the RVs from the dat file for correction:
    filename = targets{t};
    M = load(filename);
    
    % skip empty files:
    if isempty(M)
        warning([starname '.vel is an empty file']);
        continue
    end
    
    % read matrix into vectors:
    M = sortrows(M,1);
    bjd_tmp = M(:,1)-0.5;
    M(bjd_tmp<jd_min,:)=[]; % Remove measurements prior to jd_min:
    bjd = M(:,1)-0.5;
    modbjd_t = mod(bjd,1)+ObsCoo;
    rvc = M(:,2);
    rvc_err = M(:,3);
    col4 = M(:,4);
    col5 = M(:,5);
    col6 = M(:,6);
    col7 = M(:,7);
    corr_t = nan(size(bjd));
    corr_err = nan(size(bjd));
    
    %==========================================================
    % subtract the nightly mean and co-add its error:
    nights = floor(bjd);
    columns = nights-jd_min+1;
    for i=1:length(bjd)
        if drift_corr
            corr_t(i) = nights_mean(columns(i)) + (b1+modbjd_t(i)*a1);
        else
            corr_t(i) = nights_mean(columns(i));
        end
        rvc(i) = rvc(i) - corr_t(i);
        if drift_corr
            corr_err(i) = sqrt(nights_mean_err(columns(i))^2 + db1^2 + (modbjd_t(i)*da1)^2);
        else
            corr_err(i) = nights_mean_err(columns(i));
        end
        rvc_err(i) = sqrt(rvc_err(i)^2 + corr_err(i)^2);
    end
    %==========================================================
    
    % calculate the std of corrected RVs:
    bias = 1-1/4/length(rvc);
    rv_std = 1.48*mad(rvc,1)/bias;
    rv_std_vec_corr(t) = rv_std;
    rv_med_err_corr(t) = nanmedian(rvc_err);
    Nrv(t) = length(rvc);

   if savedata
    % Write an entry to .orb file (corrected RVs, NO outliers removed):
    fprintf(fid,'STAR: %s\n',starname);
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid,'A %f %f %f\n',bjd(i)+0.5,0.001*rvc(i),0.001*rvc_err(i));
        end
    end
    fprintf(fid,'END\n');
    
    % Write an entry to all_rvs file (corrected RVs, NO outliers removed):
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid_all,'%14s\t%13.5f\t%9.3f\t%7.3f\t%7.4f\t%6.4f\t%6d\t%4d\t%5.3f\t%5.3f\n',starname,bjd(i)+0.5,rvc(i),rvc_err(i),col4(i),col5(i),col6(i),col7(i),corr_t(i),corr_err(i));
        end
    end
    
    % Write a .avc.vel file (corrected RVs, NO outliers removed):
    fid_avc = fopen([dirname 'corrected_rvs/' starname '.avc.vel'],'w');
    for i=1:length(bjd)
        if ~isnan(rvc(i))
            fprintf(fid_avc,'%f %f %f\n',bjd(i)+0.5,rvc(i),rvc_err(i));
        else
            Nnoaverage(t) = Nnoaverage(t)+1;
        end
    end
    fclose(fid_avc);
   end

end % end of loop on all targets

if savedata
    fclose(fid);
    fclose(fid_all);
end

end % end of unselfbias

%% ===========================================================================
% STARTING POINT OF SAVING THE RESULTS
%===========================================================================

% re-detect RV-loud stars for writing purposes
ind = find(rv_std_vec_corr>final_rv_var);

if savedata
% all stars
    fid_txt = fopen([dirname 'txt_files\hires_RVstd_' day_str '.txt'],'w');
    fprintf(fid_txt,'starname Nrv RVstd_bef RVmederr_bef RVstd_aft RVmederr_aft\n');
    for i=1:length(Nrv)
        fprintf(fid_txt,'%s\t%d\t%7.2f\t%7.2f\t%7.2f\t%7.2f\n',stars{i},Nrv(i),rv_std_vec(i),rv_med_err(i),rv_std_vec_corr(i),rv_med_err_corr(i));
    end
    fclose(fid_txt);

% make a txt table for RV-loud stars
    fid_act = fopen([dirname 'txt_files\rv-loud_stars_' day_str '.txt'],'w');
    fprintf(fid_act,'starname   Nrv    RVstd RVstd_err RV_mederr\n');
    for i=1:length(ind)
        fprintf(fid_act,'%s\t%d\t%7.1f\t%7.1f\t%7.1f\n',stars{ind(i)},Nrv(ind(i)),rv_std_vec_corr(ind(i)),rv_std_vec_corr(ind(i))/sqrt(2*Nrv(ind(i))),rv_med_err_corr(ind(i)));
    end
    fclose(fid_act);

% make a txt table for the Nightly averages
    fid_nav = fopen([dirname 'txt_files\hires_night_zero_' day_str '.txt'],'w');
    fprintf(fid_nav,'night Nrv zero err\n');
    for i=1:length(nights_mean)
        fprintf(fid_nav,'%d\t%d\t%5.2f\t%4.2f\n',jd_min+nights_vec(i)-1,Nrv_night(i),nights_mean(i),nights_mean_err(i));
    end
    fclose(fid_nav);

% save the output data
    save([dirname 'mat_files\hires_stars_' day_str],'stars','Nrv','rvc_mean_vec','rvc_mean_vec_err','Nnoaverage','Nout_star','Nout_night','ind_same_vec','rv_std_vec','rv_med_err','nights_mean','nights_mean_err','nights_mean_std','Nrv_night','rv_std_vec_corr','rv_med_err_corr','jd_mat','rv_mat','err_mat','nzp','enzp');
end

%%
toc