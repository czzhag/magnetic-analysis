t1={'220114 02:14:06','220114 03:16:17','220114 04:16:55','220114 05:17:32','220114 06:17:34','220114 07:18:11','220114 08:18:48'};
t2={'220114 03:16:15','220114 04:16:54','220114 05:17:31','220114 06:17:33','220114 07:18:10','220114 08:18:47','220114 09:27:48'};
mydir='~/analysis/ba/ba/misc/cheng_analysis/202201_magBA/mag_coil_bias';

% process the coil data and merge results into one data file
if ~exist(sprintf('%s/magcoil_dat.mat',mydir))

  A1=[];A2=[];PHI1=[];PHI2=[];bias=[];lc=struct();lcdat=[];
  for i=1:length(t1)
  
    % process data if result doesn't exist
    if ~exist(sprintf('%s/magcoil_dat_%d.mat',mydir,i))
      analyze_mag_bias(t1{i},t2{i},mydir,sprintf('%d',i));
    end
  
    % load result and merge
    disp(sprintf('Load exist data: %s/magcoil_dat_%d.mat',mydir,i));
    d=load(sprintf('%s/magcoil_dat_%d.mat',mydir,i));
  
    A1=cat(1,A1,d.A1);
    A2=cat(1,A2,d.A2);
    PHI1=cat(1,PHI1,d.PHI1);
    PHI2=cat(1,PHI2,d.PHI2);
    bias=cat(2,bias,d.bias);
    
    % find the real load curves
    lcflns=fieldnames(d.lc);
    rglc=find(any(~isnan(d.lc.chisq),2));
   
    if ~isempty(rglc)
      for i=1:length(lcflns)
        fln=lcflns{i};
        if size(d.lc.(fln),1)==length(d.lc.s)
          if ~isfield(lc,fln)
            lc.(fln)=[];
          end
          lc.(fln)=cat(1,lc.(fln),d.lc.(fln)(rglc,:));
        else
          lc.(fln)=d.lc.(fln);
        end
      end
  
      lcdat=cat(1,lcdat,d.lcdat(rglc,:));
    end
      
  end 
  
  fn=sprintf('%s/magcoil_dat.mat',mydir);
  save(fn,'A1','A2','PHI1','PHI2','bias','lc','lcdat','-v7.3');

end

% load the merged results
load(sprintf('%s/magcoil_dat.mat',mydir))
[p ind]=get_array_info('20220101');
[fbu_calfac bias_calfac]=get_calfac(p);

% !! make some plots !!
figdir=sprintf('%s/plots',mydir);
if ~exist(figdir,'dir')
  mkdir(figdir)
end

% col 10,11 and col 16,17 are AUX channels in the first and second half of the schedule. The purpose was to monitor the phase. That's something we can worry about later. Now I will mask out the AUX channels, and use the nanmean of the amplitudes from the two halves.
% This part involves hard coding, like the 39 bias points is specifically for these runs took at the beginning of 2022.
A1m=nan(2,39,size(A1,2));A1m(1,:,:)=A1(1:39,:);A1m(2,:,:)=A1(40:78,:);
A1m(1,:,find(p.phys_rx==0&ismember(p.mce_col,[10 11])))=nan;
A1m(2,:,find(p.phys_rx==0&ismember(p.mce_col,[16 17])))=nan;
A1m=squeeze(nanmean(A1m,1));
A2m=nan(2,39,size(A2,2));A2m(1,:,:)=A2(1:39,:);A2m(2,:,:)=A2(40:78,:);
A2m(1,:,find(p.phys_rx==0&ismember(p.mce_col,[10 11])))=nan;
A2m(2,:,find(p.phys_rx==0&ismember(p.mce_col,[16 17])))=nan;
A2m=squeeze(nanmean(A2m,1));
bl=bias(1:39);

% you can first check if the coil measurement is consistant with the plt on elnod/azs
% The two schedules are run close in both time and configuration (
% 1_ba1_magnoise_1_001.sch 220114 02:14:06 - 220114 09:27:48  
% 1_ba_magnoise_responsivity_1_001.sch 220112 22:25:00 - 220113 03:40:16).
% I will assume they have very close loading for now and may check the load curve later. 
noisedir = '~/analysis/ba/ba/misc/cheng_analysis/202201_magBA';
d_plt1 = load(sprintf('%s/magnoise_plt1_lcdat.mat',noisedir));
g_plt1 = load(sprintf('%s/magnoise_plt1.mat',noisedir));
d_sky1 = load(sprintf('%s/magnoise_sky1_lcdat.mat',noisedir));
g_sky1 = load(sprintf('%s/magnoise_sky1.mat',noisedir));
n_sky1 = load(sprintf('%s/net_bias_sky1/net_20220117/noise_analysis.mat',noisedir));
d_sky2 = load(sprintf('%s/magnoise_sky2_lcdat.mat',noisedir));
g_sky2 = load(sprintf('%s/magnoise_sky2.mat',noisedir));
indplot = find(ismember(1:792,ind.rgl));

% First look of coil/mag elnod/mag az response
% All three were measured with the blanking plate on. (So loading should be similar).
% Make 3 mesh plots. X shows different bias (39 in total, from high to low).
% Plot saved in analysis_data/202201_magBA/mag_coil_bias/plots/coil_vs_pltonenazs.png
if 0
fig=figure('position',[10 10 1200 1000]);
subplot(1,3,1)
absA1m=abs(A1m);
normabsA1m=absA1m./repmat(nanmax(absA1m,[],1),size(absA1m,1),1);
imagesc(normabsA1m(:,indplot)');caxis([0 1]);
title('BA1 coil response','fontsize',14);
set(gca,'fontsize',14);
subplot(1,3,2)
abseng=abs(g_plt1.en.g(:,:,2));
normabseng=abseng./repmat(nanmax(abseng,[],1),size(abseng,1),1);
imagesc(normabseng(:,indplot)');caxis([0 1]);
title('BA1 elnod with blanking plate','fontsize',14);
set(gca,'fontsize',14);
subplot(1,3,3)
absazs=abs(g_plt1.azs.c(:,:,2));
normabsazs=absazs./repmat(nanmax(absazs,[],1),size(absazs,1),1);
imagesc(normabsazs(:,indplot)');caxis([0 1]);
title('BA1 azscan with blanking plate','fontsize',14);
set(gca,'fontsize',14);

fn=sprintf('%s/coil_vs_pltonenazs.png',figdir);
print(fig,fn,'-dpng')
disp(sprintf('made plot: %s',fn));
end

% check response
% Plots saved in analysis_data/202201_magBA/mag_coil_bias/plots/respb_check/en_nl_gcp*.png
% Left: peak normalized NEI, on-sky elnod, coil response. Mag elnod is normalized with the peak value of the on-sky elnod,
%	vs 1/v, approximately the responsivity. 
% Right: Scatter of NEI vs on-sky elnod (slope is something similar to NET)
% Conclusions: #TODO
if 0

figdirsb = sprintf('%s/respb_check',figdir);
if ~exist(figdirsb,'dir')
  mkdir(figdirsb);
end

for i = 1:length(indplot)

  j=indplot(i);
  if isempty(d_sky1.lcdat(j).rdet)
    continue
  end
  bl_cal=bl*bias_calfac(j);

  rl_sky1 = interp1(d_sky1.lcdat(j).bias(end:-1:1),d_sky1.lcdat(j).rdet(end:-1:1),bl_cal);
  inv_sky1=(p.r_sh(j)+rl_sky1)./(p.r_sh(j).*rl_sky1)./bl_cal;
  ptrans = find(rl_sky1>0.02*rl_sky1(1)&rl_sky1<0.98*rl_sky1(1));

  rl_plt1 = interp1(d_plt1.lcdat(j).bias(end:-1:1),d_plt1.lcdat(j).rdet(end:-1:1),bl_cal);
  inv_plt1=(p.r_sh(j)+rl_plt1)./(p.r_sh(j).*rl_plt1)./bl_cal;
  ptrans_plt1 = find(rl_plt1>0.02*rl_plt1(1)&rl_plt1<0.98*rl_plt1(1));

  if isempty(ptrans) | isempty(ptrans_plt1)
    continue
  end

  fig=figure('position',[10 10 1200 600]);

  subplot(1,2,1)
  plot(inv_sky1(ptrans),g_sky1.en.g(ptrans,j,2)/nanmax(g_sky1.en.g(ptrans,j,2)),'b');hold all;
  plot(inv_plt1(ptrans_plt1),g_plt1.en.g(ptrans_plt1,j,2)/nanmax(g_sky1.en.g(ptrans,j,2)),'r');
  plot(inv_sky1(ptrans),n_sky1.nl_uncal(ptrans,j)/nanmax(n_sky1.nl_uncal(ptrans,j)),'k');
  plot(inv_plt1(ptrans_plt1),A1m(ptrans_plt1,j)./nanmax(A1m(ptrans_plt1,j)),'g');

  x=1:length(inv_sky1(ptrans));y=inv_sky1(ptrans);
  c=y/nanmax(y)*10;
  scatter(inv_sky1(ptrans),g_sky1.en.g(ptrans,j,2)/nanmax(g_sky1.en.g(ptrans,j,2)),[],c,'marker','+');
  %scatter(inv_plt1(ptrans_plt1),g_plt1.en.g(ptrans_plt1,j,2)/nanmax(g_sky1.en.g(ptrans,j,2)),[],c,'marker','d');
  scatter(inv_sky1(ptrans),n_sky1.nl_uncal(ptrans,j)/nanmax(n_sky1.nl_uncal(ptrans,j)),[],c,'marker','.');
  %scatter(inv_plt1(ptrans_plt1),A1m(ptrans_plt1,j)./nanmax(A1m(ptrans_plt1,j)),60,'g');
  colormap jet  

  legend({'normed elnod gain', 'normed magnetic elnod (against sky elnod)', 'normed current noise','normed coil response'},'fontsize',14);
  
  xlim([0 0.2e8]);ylim([-0.2 1.2]);
  xlabel('1/V', 'fontsize',14); 
  ylabel('normed elnod/noise/coil', 'fontsize',14);
  set(gca,'fontsize',14);
  title(sprintf('gcp %d, tile %d',p.gcp(j),p.tile(j)),'fontsize',14);
  grid on;
  box on;

  subplot(1,2,2)
  scatter(g_sky1.en.g(ptrans,j,2)/nanmax(g_sky1.en.g(ptrans,j,2)),n_sky1.nl_uncal(ptrans,j)/nanmax(n_sky1.nl_uncal(ptrans,j)),[],c);
  colormap jet  
  
  xlim([-0.2 1.2]);ylim([-0.2 1.2]);
  xlabel('normed elnod','fontsize',14);
  ylabel('normed current noise','fontsize',14);
  grid on; box on;

  suptitle(sprintf('elnod, current noise, gcp%d', p.gcp(j)));
  fn=sprintf('%s/en_nl_gcp%d.png',figdirsb,p.gcp(j));
  print(fig,fn,'-dpng');

end

end

% Make plots to show NET and magnetic response
% Doing double y axes: left for NET, right for coil(or mag-az)/NEI response.
% datasets:
% NET      - sky1/2 (n_sky1.nl(39x2376))
% coil     - coil data set (plate on) A1m (39x792)
% NEI      - sky1/2 (n_sky.nl_uncal(39x2376))
% AZ       - plt1 (g_plt1.azs.c(:,:,2)(39x2376))
% Different datasets are aligned with R in X direction (Here we assume the magnetic response only dependence on R)
% lcdat:
% sky1/2   - d_sky1/2.lcdat(2376).xxx
% plt1     - d_plt1.lcdat(2376).xxx
% I'm making this plot for each detector (in rgl) for now
if 1

figdirsb = sprintf('%s/net_mag_r');
if ~exist(figdirsb,'dir')
  mkdir(figdirsb);
end



for i=1:length(indplot)
  
  j=indplot(i);

end

end
