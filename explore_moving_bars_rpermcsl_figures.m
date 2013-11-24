% explore_moving_bars
clc;
clear all;
close all;

database_select_samples; % get the right stuff...
pdatapath = 'E:\damon\new_TP_pdata';
addpath(pdatapath);

f_mb_L2_562 = find(i_L2.*i_mb_rpermcsl.*~i_moving.*i_562.*i_21Dhh); 

fL2_mb_562 = create_neuron_structure_all(f_mb_L2_562);
clear fname; clear datetime; clear stimcode; clear activecells; clear driver;

n_series = length(fL2_mb_562);
part1 = fL2_mb_562(1:round(n_series/2));
part2 = fL2_mb_562((round(n_series/2)+1):end);
clear fL2_mb_562;
save('movebar_rpermcsl_fL2_parts','part1','part2');

clear part2;
pack;
L2_mb_562_p1 = load_neuron_data5(part1,pdatapath);
save('movebar_rpermcsl_L2_p1','L2_mb_562_p1');

clear part1;
clear L2_mb_562_p1;
load('movebar_rpermcsl_fL2_parts','part2');
L2_mb_562_p2 = load_neuron_data5(part2,pdatapath);
save('movebar_rpermcsl_L2_p2','L2_mb_562_p2');
clear L2_mb_562_p2; clear part2; pack;

load('movebar_rpermcsl_L2_p1'); % L2_mb_562_p1
new_n = length(L2_mb_562_p1);
L2_mb_562_p11 = L2_mb_562_p1(1:round(new_n/4));
L2_mb_562_p12 = L2_mb_562_p1((round(new_n/4)+1):round(new_n/2));
L2_mb_562_p13 = L2_mb_562_p1((round(new_n/2)+1):(round(new_n/4)*3));
L2_mb_562_p14 = L2_mb_562_p1(((round(new_n/4)*3)+1):end);
save('movebar_rpermcsl_L2_p1_parts','L2_mb_562_p11','L2_mb_562_p12','L2_mb_562_p13','L2_mb_562_p14');
clear L2_mb_562_p1; clear L2_mb_562_p11; clear L2_mb_562_p12; clear L2_mb_562_p13; clear L2_mb_562_p14; pack;

%%

load('movebar_rpermcsl_L2_p1_parts','L2_mb_562_p11');
aL2_562_p11 = aggregate_mb_rpermcsl(L2_mb_562_p11);
save('movebar_rpermcsl_aL2_p11','aL2_562_p11');

clear L2_mb_562_p11; clear aL2_562_p11; pack;

load('movebar_rpermcsl_L2_p1_parts','L2_mb_562_p12');
aL2_562_p12 = aggregate_mb_rpermcsl(L2_mb_562_p12);
save('movebar_rpermcsl_aL2_p12','aL2_562_p12');

clear L2_mb_562_p11; clear aL2_562_p11; pack;

%%

ang = [0, pi, pi/2, 3*pi/2];
angs = ang;
[sortang,sortord]=sort(ang);
sigangcomp = [];

%% Select traces that worked - 562

real_inds = zeros(length(L2_mb_562),1);
inverted_inds = zeros(length(L2_mb_562),1);
for ii = 1:length(L2_mb_562)
    
    s = L2_mb_562(ii).stim; 
    ch = find(diff(s));
    ch = [1; ch(1:7)]; % keep first 7, add the 1st index...
    sig = L2_mb_562(ii).ratio;
    norm = mean(sig);
    sig = sig/norm;
    clear tr;
    % kluge for movement artifacts
    for i=1:8
        if((ch(i)+min(diff(ch)))<=length(sig))
            trtemp=sig(:,ch(i):ch(i)+min(diff(ch)));
        else
            trtemp=sig(:,ch(i):end);
        end
        trtemp(trtemp>1.6)=median(trtemp(:));
        trtemp(trtemp<0.4)=median(trtemp(:));
        tr{i}=trtemp;
    end
    t_tr=L2_mb_562(ii).t(1:size(tr{i},2)); %/in.xml.framerate;

    clear tr_n mpt xc yc rc sigang;
    tr_n = zeros(8,length(t_tr));
    for j=1:8
        tr_n(j,:)=tr{j}-mean(tr{j});
    end

    % t_tr = time trace for one pass of the bar across the screen
    % tr = cell array, length = 8 (number of different angles)
    % tr{i} = matrix of size n_cells * n_tpts
    % tr_n = same as tr, normalized by mean, mean subtracted

    [B,A]=butter(1,1/4,'low');
    % method of maximum cross correlation to find the center of the RF for
    % each cell
    % try filtering first over 5?
    trtemp=filtfilt(B,A,tr_n');
    [c,lags]=xcorr(trtemp,'coeff');
    [rc(1),mpt(1)]=max(c(:,2)); % 1 vs. 2
    [rc(2),mpt(2)]=max(c(:,2*8+4)); % 3 vs. 4
    [rc(3),mpt(3)]=max(c(:,4*8+6)); % 5 vs. 6
    [rc(4),mpt(4)]=max(c(:,6*8+8)); % 7 vs. 8

    mpt=lags(mpt); % frames off...
    posmax=mpt/length(t_tr)/2*90*sqrt(2);
    meanrc=mean(rc,2);

    % now, solve for best x,y fit to equations
    % eq 1: x = posmax(1)
    % eq 2: y = posmax(2)
    % eq 3: (x+y)/sqrt(2) = posmax(3)
    % eq 4: (x-y)/sqrt(2) = posmax(4)
    A = [1 0; ...
        0 -1;  ...
        2^-.5 2^-.5; ...
        2^-.5 -2^-.5];
    b = posmax(:);
    x0=(A'*A)^-1*A'*b;
    xc=x0(1);
    yc=x0(2);
    locpos=round(A*x0/(sqrt(2)*90)*length(t_tr)+length(t_tr)/2); % time position of max
    locpos8=[locpos(1),length(t_tr)-locpos(1),...
        locpos(2),length(t_tr)-locpos(2),...
        locpos(3),length(t_tr)-locpos(3),...
        locpos(4),length(t_tr)-locpos(4)];

    [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
    x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);

    % Show results
    figure; hold on;
    for i=1:8
        R=repmat(tr{i},[length(t_tr),1]);
        x1=x*cos(angs(i))+y*sin(angs(i));
        y1=-x*sin(angs(i))+y*cos(angs(i));
        mesh(x1,y1,R);
    end
    xlabel('x pos');
    ylabel('y pos');
    zlabel('dR/R');

    plot3(xc,yc,2,'k*','markersize',12);

    if(strcmp(input('real?','s'),'y'))
        real_inds(ii) = 1;
        if(strcmp(input('inverted?','s'),'y'))
            inverted_inds(ii) = 1;
        end
        if(~strcmp(input('correct center?','s'),'y'))

            %         Manual selection of cell centers
            figure; hold on;
            for i=1:8
                R=repmat(tr{i},[length(t_tr),1]);
                x1=x*cos(angs(i))+y*sin(angs(i));
                y1=-x*sin(angs(i))+y*cos(angs(i));
                mesh(x1,y1,R);
            end
            xlabel('x pos');
            ylabel('y pos');
            zlabel('dR/R');

            [xc,yc]=ginput(1);
            plot3(xc,yc,2,'wo');
            plot3(xc,yc,2,'w*');
            x0 = [xc; yc];

            locpos=round(A*x0/(sqrt(2)*90)*length(t_tr)+length(t_tr)/2); % time position of max
            locpos8=[locpos(1),length(t_tr)-locpos(1),...
                locpos(2),length(t_tr)-locpos(2),...
                locpos(3),length(t_tr)-locpos(3),...
                locpos(4),length(t_tr)-locpos(4)];
        end
    end

    L2_mb_562(ii).locations = locpos8;
    L2_mb_562(ii).center = x0;

    %     drawnow;
    %     pause;
    close all;
    %     clear mtr;
end
real_inds_562 = real_inds;
inverted_inds_562 = inverted_inds;
save('L2_562_inds','real_inds_562','inverted_inds_562','bad_562_inds','L2_mb_562');

%%

load('L2_562_inds');
% load('L2_562_select_inds');

aL2ni_562 = aL2_562;
aL2ni_562.neuron = aL2ni_562.neuron((real_inds_562==1)&(inverted_inds_562==0)&(bad_562_inds==0));
aL2ni_562.flyID = aL2ni_562.flyID((real_inds_562==1)&(inverted_inds_562==0)&(bad_562_inds==0));
aL2ni_562.name = aL2ni_562.name((real_inds_562==1)&(inverted_inds_562==0)&(bad_562_inds==0));
for k = 1:8
    aL2ni_562.rats{k} = aL2ni_562.rats{k}((real_inds_562==1)&(inverted_inds_562==0)&(bad_562_inds==0),:);
end

% aL2ni_562_select = aL2_562_select;
% aL2ni_562_select.neuron = aL2ni_562_select.neuron((real_inds_562_select==1)&(inverted_inds_562_select==0)&(bad_562_select_inds==0));
% aL2ni_562_select.flyID = aL2ni_562_select.flyID((real_inds_562_select==1)&(inverted_inds_562_select==0)&(bad_562_select_inds==0));
% aL2ni_562_select.name = aL2ni_562_select.name((real_inds_562_select==1)&(inverted_inds_562_select==0)&(bad_562_select_inds==0));
% for k = 1:8
%     aL2ni_562_select.rats{k} = aL2ni_562_select.rats{k}((real_inds_562_select==1)&(inverted_inds_562_select==0)&(bad_562_select_inds==0),:);
% end

L2ni_mb_562 = L2_mb_562((real_inds_562==1)&(inverted_inds_562==0)&(bad_562_inds==0));
L2_ni_inds_562 = find((real_inds_562==1)&(inverted_inds_562==0)&(bad_562_inds==0));
% L2ni_mb_562_select = L2_mb_562_select((real_inds_562_select==1)&(inverted_inds_562_select==0)&(bad_562_select_inds==0));
% L2_ni_inds_562_select = find((real_inds_562_select==1)&(inverted_inds_562_select==0)&(bad_562_select_inds==0));

aL2i_562 = aL2_562;
aL2i_562.neuron = aL2i_562.neuron((real_inds_562==1)&(inverted_inds_562==1)&(bad_562_inds==0));
aL2i_562.flyID = aL2i_562.flyID((real_inds_562==1)&(inverted_inds_562==1)&(bad_562_inds==0));
aL2i_562.name = aL2i_562.name((real_inds_562==1)&(inverted_inds_562==1)&(bad_562_inds==0));
for k = 1:8
    aL2i_562.rats{k} = aL2i_562.rats{k}((real_inds_562==1)&(inverted_inds_562==1)&(bad_562_inds==0),:);
end

% aL2i_562_select = aL2_562_select;
% aL2i_562_select.neuron = aL2i_562.neuron((real_inds_562_select==1)&(inverted_inds_562_select==1)&(bad_562_select_inds==0));
% aL2i_562_select.flyID = aL2i_562.flyID((real_inds_562_select==1)&(inverted_inds_562_select==1)&(bad_562_select_inds==0));
% aL2i_562_select.name = aL2i_562.name((real_inds_562_select==1)&(inverted_inds_562_select==1)&(bad_562_select_inds==0));
% for k = 1:8
%     aL2i_562_select.rats{k} = aL2i_562_select.rats{k}((real_inds_562_select==1)&(inverted_inds_562_select==1)&(bad_562_select_inds==0),:);
% end

L2i_mb_562 = L2_mb_562((real_inds_562==1)&(inverted_inds_562==1)&(bad_562_inds==0));
% L2i_mb_562_select = L2_mb_562((real_inds_562_select==1)&(inverted_inds_562_select==1)&(bad_562_select_inds==0));

%% Mean response as a function of angle

t = [1:size(aL2ni_562(1).rats{1},2)]/100;
mod_t = ((-t(end)/2)+(10^-2)):(10^-2):(t(end)/2);
% range = 275:1275;
% range = 425:675;
range = 1:length(t);

figure; hold on;
patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 2 2]',[0.85 0.85 0.85],...
    'EdgeColor','none');
for ii=1:8
    x = aL2ni_562.rats{sortord(ii)};
    [y dum] = align_mb2(x,L2ni_mb_562,sortord(ii));
    [m,s] = mean_cat(y,1,aL2ni_562.flyID);
    mL2ni_562(ii,:) = m; sL2ni_562(ii,:) = s;
    m = m-1; 
%     len = round(length(m)/2);
%     m = [m((len+1):end) m(1:len)];
%     s = [s((len+1):end) s(1:len)];
    plot_err_patch(mod_t(range),(ii-1)/4+m(range),s(range),[0 0 1],[.5 .5 1]);
    line([mod_t(range(1)) mod_t(range(end))],[(ii-1)/4 (ii-1)/4],...
        'color',[0 0 0],'LineStyle','--');
    text(0.5+mod_t(range(1)),0.05+(ii-1)/4,num2str(sortang(ii)*180/pi));
    title(['L2 moving bar responses, 562nm, aligned (N flies = '...
        num2str(length(unique(aL2ni_562.flyID))) ', N cells = ' num2str(length(aL2ni_562.neuron)) ')']);
end
xlabel('time (sec)');
ylabel('response (dR/R)');
ylim([-0.15 2]);xlim([mod_t(range(1)) mod_t(range(end))]);
single = [-0.1 0 0.1];
ticks = zeros(8*3,1);
for i = 1:8
    ticks(((i-1)*3+1):i*3) = single+(i-1)/4;
end
set(gca,'yTick',ticks);
set(gca,'yTickLabel',...
    '-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|');
niceaxes;

% figure; hold on;
% patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 2 2]',[0.85 0.85 0.85],...
%     'EdgeColor','none');
% for ii=1:8
%     x = aL2ni_562_select.rats{sortord(ii)};
%     [y dum] = align_mb2(x,L2ni_mb_562_select,sortord(ii));
%     [m,s] = mean_cat(y,1,aL2ni_562_select.flyID);
%     mL2ni_562_select(ii,:) = m; sL2ni_562_select(ii,:) = s;
%     m = m-1; 
%     len = round(length(m)/2);
%     m = [m((len+1):end) m(1:len)];
%     s = [s((len+1):end) s(1:len)];
%     plot_err_patch(mod_t(range),(ii-1)/4+m(range),s(range),[0 0 1],[.5 .5 1]);
%     line([mod_t(range(1)) mod_t(range(end))],[(ii-1)/4 (ii-1)/4],...
%         'color',[0 0 0],'LineStyle','--');
%     text(0.5+mod_t(range(1)),0.05+(ii-1)/4,num2str(sortang(ii)*180/pi));
%     title(['L2 moving bar responses, 562_selectnm, aligned (N flies = '...
%         num2str(length(unique(aL2ni_562_select.flyID))) ', N cells = ' num2str(length(aL2ni_562_select.neuron)) ')']);
% end
% xlabel('time (sec)');
% ylabel('response (dR/R)');
% ylim([-0.15 2]);xlim([mod_t(range(1)) mod_t(range(end))]);
% single = [-0.1 0 0.1];
% ticks = zeros(8*3,1);
% for i = 1:8
%     ticks(((i-1)*3+1):i*3) = single+(i-1)/4;
% end
% set(gca,'yTick',ticks);
% set(gca,'yTickLabel',...
%     '-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|-0.1|0|0.1|');
% niceaxes;

%% Mean response as a function of angle - without ticks

t = [1:size(aL2ni_562(1).rats{1},2)]/100;
mod_t = ((-t(end)/2)+(10^-2)):(10^-2):(t(end)/2);
range = 275:1275;

figure; hold on;
patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 2 2]',[0.85 0.85 0.85],...
    'EdgeColor','none');
mL2ni_562 = zeros(8,size(aL2ni_562.rats{1},2));
sL2ni_562 = zeros(8,size(aL2ni_562.rats{1},2));
for ii=1:8
    x = aL2ni_562.rats{sortord(ii)};
    [y dum] = align_mb2(x,L2ni_mb_562,sortord(ii));
    [m,s] = mean_cat(y,1,aL2ni_562.flyID);
    m = m-1; 
%     len = round(length(m)/2);
%     m = [m((len+1):end) m(1:len)];
%     s = [s((len+1):end) s(1:len)];
    mL2ni_562(ii,:) = m; sL2ni_562(ii,:) = s;
    plot_err_patch(mod_t(range),(ii-1)/4+m(range),s(range),[0 0 1],[.5 .5 1]);
    line([mod_t(range(1)) mod_t(range(end))],[(ii-1)/4 (ii-1)/4],...
        'color',[0 0 0],'LineStyle','--');
    text(0.5+mod_t(range(1)),0.07+(ii-1)/4,num2str(sortang(ii)*180/pi));
    title(['L2 moving bar responses, 562nm, aligned (N flies = '...
        num2str(length(unique(aL2ni_562.flyID))) ', N cells = ' num2str(length(aL2ni_562.neuron)) ')']);
end
xlabel('time (sec)');
ylabel('response (dR/R)');
ylim([-0.15 2]);xlim([mod_t(range(1)) mod_t(range(end))]);
line([1.5+mod_t(range(1)) 1.5+mod_t(range(1))],[1.8 1.9],'color',[0 0 0]);
line([1.5+mod_t(range(1)) 2+mod_t(range(1))],[1.8 1.8],'color',[0 0 0]);
niceaxes;

% figure; hold on;
% patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 2 2]',[0.85 0.85 0.85],...
%     'EdgeColor','none');
% mL2ni_562_select = zeros(8,size(aL2ni_562_select.rats{1},2));
% sL2ni_562_select = zeros(8,size(aL2ni_562_select.rats{1},2));
% for ii=1:8
%     x = aL2ni_562_select.rats{sortord(ii)};
%     [y dum] = align_mb2(x,L2ni_mb_562_select,sortord(ii));
%     [m,s] = mean_cat(y,1,aL2ni_562_select.flyID);
%     m = m-1; 
%     len = round(length(m)/2);
%     m = [m((len+1):end) m(1:len)];
%     s = [s((len+1):end) s(1:len)];
%     mL2ni_562_select(ii,:) = m; sL2ni_562_select(ii,:) = s;
%     plot_err_patch(mod_t(range),(ii-1)/4+m(range),s(range),[0 0 1],[.5 .5 1]);
%     line([mod_t(range(1)) mod_t(range(end))],[(ii-1)/4 (ii-1)/4],...
%         'color',[0 0 0],'LineStyle','--');
%     text(0.5+mod_t(range(1)),0.07+(ii-1)/4,num2str(sortang(ii)*180/pi));
%     title(['L2 moving bar responses, 562nm select, aligned (N flies = '...
%         num2str(length(unique(aL2ni_562_select.flyID))) ', N cells = ' num2str(length(aL2ni_562_select.neuron)) ')']);
% end
% xlabel('time (sec)');
% ylabel('response (dR/R)');
% ylim([-0.15 2]);xlim([mod_t(range(1)) mod_t(range(end))]);
% line([1.5+mod_t(range(1)) 1.5+mod_t(range(1))],[1.8 1.9],'color',[0 0 0]);
% line([1.5+mod_t(range(1)) 2+mod_t(range(1))],[1.8 1.8],'color',[0 0 0]);
% niceaxes;

%% direction selectivity?

t = [1:size(aL2ni_562(1).rats{1},2)]/100;
mod_t = ((-t(end)/2)+(10^-2)):(10^-2):(t(end)/2);
range = 275:1275;

directions = {'diag1','horizontal','diag2','vertical'};

for ii = 1:4
    figure; hold all;
    patch([-1.6 -1.4 -1.4 -1.6]',[-0.2 -0.2 0.1 0.1]',[0.85 0.85 0.85],...
        'EdgeColor','none');
    m1 = mL2ni_562(ii,:);
    s1 = sL2ni_562(ii,:);
    m2 = mL2ni_562(ii+4,:);
    s2 = sL2ni_562(ii+4,:);
    [temp,delays] = xcorr(m1,m2);
    ind = delays(find(temp==max(temp)));
    m2 = circshift(m2',ind)';
    s2 = circshift(s2',ind)';
    plot_err_patch_v2(mod_t(range),m1(range),s1(range),[0 0 1],[.5 .5 1]);
    plot_err_patch_v2(mod_t(range),m2(range),s2(range),0.5*[0 0 1],0.5*[.5 .5 1]);
    line([mod_t(range(1)) mod_t(range(end))],[0 0],'color',[0 0 0]);
    title(['L2 moving bar responses, 562nm, direction ' directions{ii}]);
    xlabel('time (sec)');
    ylabel('response (dR/R)');
    ylim([-0.2 0.1]);xlim([mod_t(range(1)) mod_t(range(end))]);
    niceaxes;
end
% for ii = 1:4
%     figure; hold all;
%     patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 0.15 0.15]',[0.85 0.85 0.85],...
%         'EdgeColor','none');
%     m1 = mL2ni_562_select(ii,:);
%     s1 = sL2ni_562_select(ii,:);
%     m2 = mL2ni_562_select(ii+4,:);
%     s2 = sL2ni_562_select(ii+4,:);
%     plot_err_patch_v2(mod_t(range),m1(range),s1(range),[0 0 1],[.5 .5 1]);
%     plot_err_patch_v2(mod_t(range),m2(range),s2(range),0.5*[0 0 1],0.5*[.5 .5 1]);
%     line([mod_t(range(1)) mod_t(range(end))],[0 0],'color',[0 0 0]);
%     title(['L2 moving bar responses, 562nm select, direction ' directions{ii}]);
%     xlabel('time (sec)');
%     ylabel('response (dR/R)');
%     ylim([-0.15 0.15]);xlim([mod_t(range(1)) mod_t(range(end))]);
%     niceaxes;
% end

%% orientation dependence ?

t = [1:size(aL2ni_562(1).rats{1},2)]/100;
mod_t = ((-t(end)/2)+(10^-2)):(10^-2):(t(end)/2);
range = 275:1275;

directions = {'diag1','horizontal','diag2','vertical'};
% directions = {'vertical','diag2','horizontal','diag1'};

figure; hold all;
patch([-1.6 -1.4 -1.4 -1.6]',[-0.2 -0.2 0.1 0.1]',[0.85 0.85 0.85],...
    'EdgeColor','none');
h = zeros(4,1);
for ii = 4:-1:1
    m1 = mL2ni_562(ii,:);
    mL2ni_562(ii+4,:);
    [temp,delays] = xcorr(m1,m2);
    ind = delays(find(temp==max(temp)));
    m2 = circshift(m2',ind)';
    s1 = sL2ni_562(ii,:);
    s2 = sL2ni_562(ii+4,:);
    s2 = circshift(s2',ind)';
    m = (m1+m2)/2;
    s = sqrt(s1.^2+s2.^2);
    h(ii) = plot_err_patch_v2(mod_t(range),m(range),s(range),0.25*ii*[0 0 1],0.25*ii*[.5 .5 1]);
end
line([mod_t(range(1)) mod_t(range(end))],[0 0],'color',[0 0 0]);
title('L2 moving bar responses, 562nm');
xlabel('time (sec)');
ylabel('response (dR/R)');
ylim([-0.2 0.1]);xlim([mod_t(range(1)) mod_t(range(end))]);
legend(h,directions,'location','southeast');
niceaxes;

% figure; hold all;
% patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 0.15 0.15]',[0.85 0.85 0.85],...
%     'EdgeColor','none');
% h = zeros(4,1);
% for ii = 1:4
%     m = (mL2ni_562_select(ii,:)+mL2ni_562_select(ii+4,:))/2;
%     s = sqrt(sL2ni_562_select(ii,:).^2+sL2ni_562_select(ii+4,:).^2);
%     h(ii) = plot_err_patch_v2(mod_t(range),m(range),s(range),0.25*ii*[0 0 1],0.25*ii*[.5 .5 1]);
% end
% line([mod_t(range(1)) mod_t(range(end))],[0 0],'color',[0 0 0]);
% title('L2 moving bar responses, 562nm select');
% xlabel('time (sec)');
% ylabel('response (dR/R)');
% ylim([-0.15 0.15]);xlim([mod_t(range(1)) mod_t(range(end))]);
% legend(h,directions);
% niceaxes;

%% mean MB response figures

range = 1:size(mL2ni_562,2);
[y,lags] = align_array_LB2(mL2ni_562,300);
m_mL2ni_562 = mean(y,1);
% m_mL2ni_562_select = mean(mL2ni_562_select,1);

x = zeros(size(sL2ni_562));
for i = 1:size(sL2ni_562,1)
    x(i,:) = circshift(sL2ni_562(i,:)',lags(i))';
end
s_sL2ni_562 = sqrt(mean(x.^2,1));
% s_sL2ni_562_select = sqrt(mean(sL2ni_562_select.^2,1));

mod_t = ((-t(end)/2)+(10^-2)):(10^-2):(t(end)/2);

figure;hold all;niceaxes;
patch([-1.6 -1.4 -1.4 -1.6]',[-0.2 -0.2 2 2]',[0.85 0.85 0.85],...
    'EdgeColor','none');
h = plot_err_patch_v2(mod_t(range),m_mL2ni_562,s_sL2ni_562,[0 0 1],[0.5 0.5 1]);
legend(h,sprintf('L2: %d (%d)',length(unique(aL2ni_562.flyID)),length(aL2ni_562.neuron)));
line([0 0],[-0.2 1.15],'color',[0 0 0],'lineWidth',1);
line([-4 7],[0 0],'color',[0 0 0],'lineWidth',1);
xlabel('time (sec)');ylabel('response (dR/R)');
ylim([-0.2 0.1]);xlim([mod_t(range(1)) mod_t(range(end))]);
title('moving bar response 562 nm');

% figure;hold all;niceaxes;
% patch([-1.6 -1.4 -1.4 -1.6]',[-0.15 -0.15 2 2]',[0.85 0.85 0.85],...
%     'EdgeColor','none');
% h = plot_err_patch_v2(mod_t(range),m_mL2ni_562_select,s_sL2ni_562_select,[1 0 0],[1 0.5 0.5]);
% legend(h,sprintf('L2: %d (%d)',length(unique(aL2ni_562_select.flyID)),length(aL2ni_562_select.neuron)));
% line([0 0],[-0.15 1.15],'color',[0 0 0],'lineWidth',1);
% line([-4 7],[0 0],'color',[0 0 0],'lineWidth',1);
% xlabel('time (sec)');ylabel('response (dR/R)');
% ylim([-0.13 0.13]);xlim([mod_t(range(1)) mod_t(range(end))]);
% title('moving bar response 562 select');

%% Save center locations and names

locs = [L2_mb_562.locations];
locs = reshape(locs',8,length(locs)/8)';

centers = [L2_mb_562.center]';

neurons = [L2_mb_562.cell];

names = cell(length(L2_mb_562),1);
for ii = 1:length(L2_mb_562)
    names{ii} = L2_mb_562(ii).name;
end

save('L2_known_locations','names','centers','locs','neurons');

%% choose examples - L2

% pdatapath='C:\Documents and Settings\Limor\My Documents\damon\pData';
% % codepath = 'C:\Documents and Settings\Limor\My Documents\damon\matlab_fly2p\Limor';
% addpath(pdatapath);
% % addpath(codepath);
% 
% ii = 1;
% % nice_examples = [];
% 
% all_distances = [];
% while ii <= length(aL2ni.neuron)
%     
%     close all;
%     sInd = ii;
%     fname = aL2ni.name{ii};
%     done = false;
%     active_cells = aL2ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL2ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL2ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL2ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1));
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L2_mb(L2_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
% %     figure; hold on;
% %     cm=colormap('lines');
% %     for i=1:Nneurons
% %         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
% %     end
% %     plot(t,stim/3);
% %     ylabel('fractional ratio change');
% %     xlabel('time (s)');
% %     title(sprintf('index =  %d',sInd));
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL2ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     mtr=zeros(8,length(t_tr));
%     
% %     for s = active_cells
% %         for i=1:2:7
% %             figure; hold on;
% %             plot(t_tr,tr{i}(s,:)','b','lineWidth',2);
% %             plot(t_tr,tr{i+1}(s,end:-1:1)','r','lineWidth',2);
% %             xlabel('+/- time (s)');
% %             ylabel('fractional change');
% %             title(['stim # ' num2str([i i+1]) ', cell = ' num2str(s)]);
% %             legend(num2str(ang(i)*180/pi),num2str(ang(i+1)*180/pi));
% %         end
% %     end
%     
%     for i=1:length(active_cells)
%         if (L2_mb(L2_ni_inds(sInd-1+i)).wavelength==490)
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1);
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2);
% %         locpos8 = L2_mb(L2_ni_inds(sInd-1+i)).locations;
%         else
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1)*2;
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2)*2+45;
%         end
%     end
%     
%     for i = 1:length(active_cells)
%         cur_center = [xc(i) yc(i)];
%         all_dists = sqrt((xc - repmat(xc(i),1,length(xc))).^2+...
%             (yc - repmat(yc(i),1,length(yc))).^2);
%         all_dists = all_dists(all_dists>0);
%         dist = min(all_dists);
%         if(dist>20)
%             disp(sprintf('error in %d',n));
%         end
%         all_distances = [all_distances; dist];
%         disp(sprintf('distance is %0.5g',dist));
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
% 
% %     figure; hold on;
% %     plot(xc,yc,'ro-');
% %     for i=1:length(xc)
% %         text(xc(i)+1,yc(i)+1,num2str(active_cells(i)));
% %         circ=2.5*exp(sqrt(-1)*[0:.01:2*pi+0.01]);
% %         plot(xc(i)+real(circ),yc(i)+imag(circ),'k:');
% %     end
% %     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
% %     xlabel('x center');
% %     ylabel('y center');
% %     title(['RF centers: ' num2str(active_cells)]);
% %     set(gca,'dataa',[1 1 1]);
% 
% %     figure;
% %     subplot(2,1,1);
% %     imagesc(in.AV1);
% %     axis image;
% %     title('image and masks');
% %     subplot(2,1,2);
% %     maskim=in.AV1;
% %     for i=1:length(active_cells)
% %         maskim(find(in.masks{active_cells(i)}))=(i+1)*((2^12)/(Nneurons+1));
% %     end
% %     hold on;imagesc(flipud(maskim));
%     
% %     drawnow;
% %     pause;
% %     close all;
%     clear mtr;
%     
% %     if(strcmp(input('nice?','s'),'y'))
% %         nice_examples = [nice_examples sInd];
% %     end
% end
% 
% % save('L2_all_distances','all_distances');
% % save('my_examples','nice_examples');
% 
% %% Making detailed figures for specific examples
% 
% load('my_examples');
% load('final_MB_curves'); 
% load('best_traces');
% load('angle_vals');
% load('angle_vals2');
% load('angle_vals3');
% 
% % n = 12: nice but RFs are too close to each other
% % n = 15: nice but few (5) cells
% % n = 18: very nice, 6 cells
% % n = 19: very nice, 6 cells (better then prev.)
% % n = 29: 8 cells, not sorted nicely in space
% % n = 34: 10 cells, some are too close together
% % n = 35: 10 cells, separation may not be sufficiently clear from traces
% % n = 37: 10 cells, separation may nclose allot be sufficiently clear from traces
% % n = 39: 10 cells, spread isn't perfect, but is good
% % n = 40: 11 cells, spread isn't perfect, but is good
% % n = 41: 11 cells, pretty nice spread
% load('my_examples');
% 
% % n > 34 for 562nm
% % 37 is an option, 38 is another option - but the spread is somewhat crazy
% % 41 if we settle on only 3 responding cells
% % 47 has 4 cells, slightly weird spread
% for n = 38
%     
%     clc;
%     n = n+1;
%     sInd = nice_examples(n);
%     ii = sInd;
%     fname = aL2ni.name{ii};
%     done = false;
%     active_cells = aL2ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL2ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL2ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL2ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1))
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L2_mb(L2_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
%     figure; hold on;
%     cm=colormap('lines');
%     for i=1:Nneurons
%         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
%     end
%     dStim = diff(stim);
%     inds = find(dStim==1);
%     text_angs = {'0','180','90','270','315','135','45','225'};
%     for i = 1:length(inds)
%         line([t(inds(i)) t(inds(i))],[0 Nneurons+1],'color',[0 0 0],'linewidth',2);
%         if(i < length(inds))
%             text(((t(inds(i))+t(inds(i+1)))/2) - 1,Nneurons+0.7,text_angs{i});
%         else
%             text((t(inds(i))+t(end)-0.2)/2 - 1,Nneurons+0.7,text_angs{i});
%         end
%     end
%     ylabel('fractional ratio change');
%     xlabel('time (s)');
%     ind = find(in.fileloc=='\');
%     title(in.fileloc((ind(end)+1):end));
%     ylim([0 Nneurons+1]);
%     niceaxes;
% 
% %     figure;plot(mean(sig,1));
%     
%     clear tr;
%     for i=1:8
%         tr{i} = aL2ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     % Potential supp. figure - how we find the RF center
% %     close all;
% %     for s = 1:length(active_cells)
% %         locpos8_temp = L2_mb(L2_ni_inds(sInd-1+s)).locations;
% %         locpos8 = zeros(8,1);
% %         for i = 1:8
% %             [val,locpos8(i)] = min(abs(t_tr-t(locpos8_temp(i))));
% %         end  
% %         total = locpos8(1)+locpos8(2);
% %         for i=1:2:7
% %             figure; hold on;
% %             plot(t_tr,tr{i}(s,:)','b','lineWidth',2);
% %             plot(t_tr,tr{i+1}(s,end:-1:1)','r','lineWidth',2);
% %             scatter(t_tr(locpos8(i)),tr{i}(s,locpos8(i)),100,[0 0 1],'filled')
% %             scatter(t_tr(total-locpos8(i+1)),tr{i+1}(s,length(t_tr)+locpos8(i+1)-total+1),...
% %                 100,[1 0 0],'filled');
% %             line([t_tr(locpos8(i)) t_tr(locpos8(i))],[0.7 1.3],'color',[0 0 0],'lineWidth',1);
% %             xlabel('+/- time (s)');
% %             ylabel('fractional change');
% %             title(['stim # ' num2str([i i+1]) ', cell = ' num2str(active_cells(s))]);
% %             legend(num2str(ang(i)*180/pi),num2str(ang(i+1)*180/pi));
% %             ylim([0.6 1.4]);
% %         end
% %     end
%     
%     mtr=zeros(8,length(t_tr));
%        
%     for i=1:length(active_cells)
%         if (L2_mb(L2_ni_inds(sInd-1+i)).wavelength==490)
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1);
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2);
% %         locpos8 = L2_mb(L2_ni_inds(sInd-1+i)).locations;
%         else
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1)*2;
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2)*2;
%         end
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
%     
%     % Create a colored map of ROIs
%     CMask = zeros(in.xml.linesperframe,in.xml.pixperline,3);
%     cm = colormap('lines');
%     nMasks = length(active_cells);
%     for i = 1:nMasks
%         curColor = cm(i,:);
%         curMask = cat(3,curColor(1).*in.masks{active_cells(i)},...
%             curColor(2).*in.masks{active_cells(i)},...
%             curColor(3).*in.masks{active_cells(i)});
%         CMask = CMask + curMask;
%     end
%     figure;
% %     subplot2(1,1,1,0,0);
% %     limit = 20:115;
% %     imshow(in.AV1(:,limit),[]);
%     imshow(in.AV1,[]);
%     axis image
%     hold on;
%     CMask = CMask./max(CMask(:));
%     h = imagesc(CMask);
% %     h = imagesc(CMask(:,limit,:));
%     set(h,'AlphaData',0.2);
% %     title('image and masks');
%     axis image
% 
%     figure; hold on;
%     for i=1:length(xc)
% %         text(xc(i)+3,yc(i),num2str(active_cells(i)));
%         circ = 2.5*exp(sqrt(-1)*[0:10^(-2):(2*pi+10^(-2))]);
%         plot(xc(i),-yc(i),'.','color',cm(i,:),'markerSize',15);
%         plot(xc(i)+real(circ),-yc(i)+imag(circ),'color',cm(i,:),'linestyle',':','linewidth',2);
%     end
%     plot(xc,-yc,'color',[0 0 0]);
%     line([-20 -20],[60 65],'color',[0 0 0],'lineWidth',2);
%     line([-20 -15],[60 60],'color',[0 0 0],'lineWidth',2);
%     xlabel('x center');ylabel('y center');
%     title('RF centers');
%     niceaxes;axis image
% %     xlim([-5 20]);ylim([0 40]); 
%     
%     clear mtr;
% 
%     done = false;
% %     close all;
%     while(~done)
%         if(strcmp(input('another direction?','s'),'y'))
%             done = false;
%         else
%             done = true;
%             break;
%         end
%         direction = input('please choose separation direction (horiz,vert,diag1,diag2)','s');
%         %         order = input('please choose filter order','s');
%         %         cutoff = input('please choose filter cutoff','s');
%         %         show_dir_separation3(t_tr,tr,direction,str2num(order),str2num(cutoff));
%         %         show_dir_separation3(t_tr,tr,direction);
%         show_dir_separation5(t_tr,tr,direction);
%     end
% 
%     %     close all;
% end
% 
% for n = 38
%     
%     clc;
%     n = n+1;
%     sInd = nice_examples(n);
%     ii = sInd;
%     fname = aL2ni.name{ii};
%     done = false;
%     active_cells = aL2ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL2ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL2ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL2ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1))
%     
%     % Only keep some of the active cells
%     active_cells = active_cells(1:6);
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L2_mb(L2_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
%     figure; hold on;
%     cm=colormap('lines');
%     for i=1:Nneurons
%         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
%     end
%     dStim = diff(stim);
%     inds = find(dStim==1);
%     text_angs = {'0','180','90','270','315','135','45','225'};
%     for i = 1:length(inds)
%         line([t(inds(i)) t(inds(i))],[0 Nneurons+1],'color',[0 0 0],'linewidth',2);
%         if(i < length(inds))
%             text(((t(inds(i))+t(inds(i+1)))/2) - 1,Nneurons+0.7,text_angs{i});
%         else
%             text((t(inds(i))+t(end)-0.2)/2 - 1,Nneurons+0.7,text_angs{i});
%         end
%     end
%     ylabel('fractional ratio change');
%     xlabel('time (s)');
%     ind = find(in.fileloc=='\');
%     title(in.fileloc((ind(end)+1):end));
%     ylim([0 Nneurons+1]);
%     niceaxes;
% 
% %     figure;plot(mean(sig,1));
%     
%     clear tr;
%     for i=1:8
%         tr{i} = aL2ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     % Potential supp. figure - how we find the RF center
% %     close all;
% %     for s = 1:length(active_cells)
% %         locpos8_temp = L2_mb(L2_ni_inds(sInd-1+s)).locations;
% %         locpos8 = zeros(8,1);
% %         for i = 1:8
% %             [val,locpos8(i)] = min(abs(t_tr-t(locpos8_temp(i))));
% %         end  
% %         total = locpos8(1)+locpos8(2);
% %         for i=1:2:7
% %             figure; hold on;
% %             plot(t_tr,tr{i}(s,:)','b','lineWidth',2);
% %             plot(t_tr,tr{i+1}(s,end:-1:1)','r','lineWidth',2);
% %             scatter(t_tr(locpos8(i)),tr{i}(s,locpos8(i)),100,[0 0 1],'filled')
% %             scatter(t_tr(total-locpos8(i+1)),tr{i+1}(s,length(t_tr)+locpos8(i+1)-total+1),...
% %                 100,[1 0 0],'filled');
% %             line([t_tr(locpos8(i)) t_tr(locpos8(i))],[0.7 1.3],'color',[0 0 0],'lineWidth',1);
% %             xlabel('+/- time (s)');
% %             ylabel('fractional change');
% %             title(['stim # ' num2str([i i+1]) ', cell = ' num2str(active_cells(s))]);
% %             legend(num2str(ang(i)*180/pi),num2str(ang(i+1)*180/pi));
% %             ylim([0.6 1.4]);
% %         end
% %     end
%     
%     mtr=zeros(8,length(t_tr));
%        
%     for i=1:length(active_cells)
%         if (L2_mb(L2_ni_inds(sInd-1+i)).wavelength==490)
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1);
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2);
% %         locpos8 = L2_mb(L2_ni_inds(sInd-1+i)).locations;
%         else
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1)*2;
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2)*2;
%         end
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
%     
%     % Create a colored map of ROIs
%     CMask = zeros(in.xml.linesperframe,in.xml.pixperline,3);
%     cm = colormap('lines');
%     nMasks = length(active_cells);
%     for i = 1:nMasks
%         curColor = cm(i,:);
%         curMask = cat(3,curColor(1).*in.masks{active_cells(i)},...
%             curColor(2).*in.masks{active_cells(i)},...
%             curColor(3).*in.masks{active_cells(i)});
%         CMask = CMask + curMask;
%     end
%     figure;
% %     subplot2(1,1,1,0,0);
%     limit = 1:115;
%     imshow(in.AV1(:,limit),[]);
% %     imshow(in.AV1,[]);
%     axis image
%     hold on;
%     CMask = CMask./max(CMask(:));
% %     h = imagesc(CMask);
%     h = imagesc(CMask(:,limit,:));
%     set(h,'AlphaData',0.2);
% %     title('image and masks');
%     axis image
% 
%     figure; hold on;
%     for i=1:length(xc)
% %         text(xc(i)+3,yc(i),num2str(active_cells(i)));
%         circ = 2.5*exp(sqrt(-1)*[0:10^(-2):(2*pi+10^(-2))]);
%         plot(xc(i),-yc(i),'.','color',cm(i,:),'markerSize',15);
%         plot(xc(i)+real(circ),-yc(i)+imag(circ),'color',cm(i,:),'linestyle',':','linewidth',2);
%     end
%     plot(xc,-yc,'color',[0 0 0]);
%     line([-20 -20],[60 65],'color',[0 0 0],'lineWidth',2);
%     line([-20 -15],[60 60],'color',[0 0 0],'lineWidth',2);
%     xlabel('x center');ylabel('y center');
%     title('RF centers');
%     niceaxes;axis image
% %     xlim([-5 20]);ylim([0 40]); 
%     
%     clear mtr;
% 
%     done = false;
% %     close all;
%     while(~done)
%         if(strcmp(input('another direction?','s'),'y'))
%             done = false;
%         else
%             done = true;
%             break;
%         end
%         direction = input('please choose separation direction (horiz,vert,diag1,diag2)','s');
%         %         order = input('please choose filter order','s');
%         %         cutoff = input('please choose filter cutoff','s');
%         %         show_dir_separation3(t_tr,tr,direction,str2num(order),str2num(cutoff));
%         %         show_dir_separation3(t_tr,tr,direction);
%         show_dir_separation5(t_tr,tr,direction);
%     end
% 
%     %     close all;
% end

% %% choose examples - L1_1
%
% pdatapath='C:\Documents and Settings\Limor\My Documents\damon\pData';
% % codepath = 'C:\Documents and Settings\Limor\My Documents\damon\matlab_fly2p\Limor';
% addpath(pdatapath);
% % addpath(codepath);
% 
% ii = 1;
% nice_examples = [];
% 
% while ii <= length(aL1_1ni.neuron)
%     
%     close all;
%     sInd = ii;
%     fname = aL1_1ni.name{ii};
%     done = false;
%     active_cells = aL1_1ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL1_1ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL1_1ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL1_1ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1));
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L1_1_mb(L1_1_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
%     figure; hold on;
%     cm=colormap('lines');
%     for i=1:Nneurons
%         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
%     end
%     plot(t,stim/3);
%     ylabel('fractional ratio change');
%     xlabel('time (s)');
%     title(sprintf('index =  %d',sInd));
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL1_1ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     mtr=zeros(8,length(t_tr));
%     
% %     for s = active_cells
% %         for i=1:2:7
% %             figure; hold on;
% %             plot(t_tr,tr{i}(s,:)','b','lineWidth',2);
% %             plot(t_tr,tr{i+1}(s,end:-1:1)','r','lineWidth',2);
% %             xlabel('+/- time (s)');
% %             ylabel('fractional change');
% %             title(['stim # ' num2str([i i+1]) ', cell = ' num2str(s)]);
% %             legend(num2str(ang(i)*180/pi),num2str(ang(i+1)*180/pi));
% %         end
% %     end
%     
%     for i=1:length(active_cells)
%         xc(i) = L1_1_mb(L1_1_ni_inds(sInd-1+i)).centers(1);
%         yc(i) = L1_1_mb(L1_1_ni_inds(sInd-1+i)).centers(2);
%         locpos8 = L1_1_mb(L1_1_ni_inds(sInd-1+i)).locations;
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
% 
%     figure; hold on;
%     plot(xc,yc,'ro-');
%     for i=1:length(xc)
%         text(xc(i)+1,yc(i)+1,num2str(active_cells(i)));
%         circ=2.5*exp(sqrt(-1)*[0:.01:2*pi+0.01]);
%         plot(xc(i)+real(circ),yc(i)+imag(circ),'k:');
%     end
%     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
%     xlabel('x center');
%     ylabel('y center');
%     title(['RF centers: ' num2str(active_cells)]);
%     set(gca,'dataa',[1 1 1]);
% 
%     figure;
%     subplot(2,1,1);
%     imagesc(in.AV1);
%     axis image;
%     title('image and masks');
%     subplot(2,1,2);
%     maskim=in.AV1;
%     for i=1:length(active_cells)
%         maskim(find(in.masks{active_cells(i)}))=(i+1)*((2^12)/(Nneurons+1));
%     end
%     hold on;imagesc(flipud(maskim));
%     
% %     drawnow;
% %     pause;
% %     close all;
%     clear mtr;
%     
%     if(strcmp(input('nice?','s'),'y'))
%         nice_examples = [nice_examples sInd];
%     end
% end
% 
% nice_examples_L1_1 = nice_examples;
% save('my_examples_L1_1','nice_examples_L1_1');
% 
% %% choose examples - L1_5
% 
% pdatapath='C:\Documents and Settings\Limor\My Documents\damon\pData';
% % codepath = 'C:\Documents and Settings\Limor\My Documents\damon\matlab_fly2p\Limor';
% addpath(pdatapath);
% % addpath(codepath);
% 
% ii = 1;
% nice_examples = [];
% 
% while ii <= length(aL1_5ni.neuron)
%     
%     close all;
%     sInd = ii;
%     fname = aL1_5ni.name{ii};
%     done = false;
%     active_cells = aL1_5ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL1_5ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL1_5ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL1_5ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1));
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L1_5_mb(L1_5_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
%     figure; hold on;
%     cm=colormap('lines');
%     for i=1:Nneurons
%         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
%     end
%     plot(t,stim/3);
%     ylabel('fractional ratio change');
%     xlabel('time (s)');
%     title(sprintf('index =  %d',sInd));
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL1_5ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     mtr=zeros(8,length(t_tr));
%     
% %     for s = active_cells
% %         for i=1:2:7
% %             figure; hold on;
% %             plot(t_tr,tr{i}(s,:)','b','lineWidth',2);
% %             plot(t_tr,tr{i+1}(s,end:-1:1)','r','lineWidth',2);
% %             xlabel('+/- time (s)');
% %             ylabel('fractional change');
% %             title(['stim # ' num2str([i i+1]) ', cell = ' num2str(s)]);
% %             legend(num2str(ang(i)*180/pi),num2str(ang(i+1)*180/pi));
% %         end
% %     end
%     
%     for i=1:length(active_cells)
%         xc(i) = L1_1_mb(L1_5_ni_inds(sInd-1+i)).centers(1);
%         yc(i) = L1_1_mb(L1_5_ni_inds(sInd-1+i)).centers(2);
%         locpos8 = L1_1_mb(L1_5_ni_inds(sInd-1+i)).locations;
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
% 
%     figure; hold on;
%     plot(xc,yc,'ro-');
%     for i=1:length(xc)
%         text(xc(i)+1,yc(i)+1,num2str(active_cells(i)));
%         circ=2.5*exp(sqrt(-1)*[0:.01:2*pi+0.01]);
%         plot(xc(i)+real(circ),yc(i)+imag(circ),'k:');
%     end
%     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
%     xlabel('x center');
%     ylabel('y center');
%     title(['RF centers: ' num2str(active_cells)]);
%     set(gca,'dataa',[1 1 1]);
% 
%     figure;
%     subplot(2,1,1);
%     imagesc(in.AV1);
%     axis image;
%     title('image and masks');
%     subplot(2,1,2);
%     maskim=in.AV1;
%     for i=1:length(active_cells)
%         maskim(find(in.masks{active_cells(i)}))=(i+1)*((2^12)/(Nneurons+1));
%     end
%     hold on;imagesc(flipud(maskim));
%     
% %     drawnow;
% %     pause;
% %     close all;
%     clear mtr;
%     
%     if(strcmp(input('nice?','s'),'y'))
%         nice_examples = [nice_examples sInd];
%     end
% end
% 
% nice_examples_L1_5 = nice_examples;
% save('my_examples_L1_5','nice_examples_L1_5');
% 
% %% Go over the nice examples only - look for angle dependence, distances
% 
% % clc;
% % clear all;
% % close all;
% 
% load('my_examples');
% load('final_MB_curves'); % 'mL2ni','mL2i','mL1_1ni','mL1_1i','mL1_5ni','mL1_5i'
% nice_examples = unique(nice_examples);
% best_examples = [];
% 
% ref = mean(mL2ni,1)-mean(mean(mL2ni));
% slow_t = [1:size(aL2ni(1).rats{1},2)]/100;
% ref_t = ((-slow_t(end)/2)+(10^-2)):(10^-2):(slow_t(end)/2);
% 
% [b,a] = butter(10,0.03);  % For smoothing
% 
% % nice_distances = [];
% % for n = 1:length(nice_examples)
% % 
% %     clc;
% %     sInd = nice_examples(n);
% %     ii = sInd;
% %     fname = aL2ni.name{ii};
% %     done = false;
% %     active_cells = aL2ni.neuron(ii);
% %     while(~done)
% %         if(~strcmp(fname,aL2ni.name{ii}))
% %             done = true;
% %         else
% %             if(ii < length(aL2ni.neuron))
% %                 ii = ii+1;
% %                 active_cells = [active_cells aL2ni.neuron(ii)];
% %             else
% %                 ii = ii-1;
% %                 done = true;
% %             end
% %         end
% %     end
% %     active_cells = active_cells(1:(end-1));
% % 
% %     in=load([pdatapath '/' fname]);
% %     in=in.strct;
% % 
% %     disp(in.fileloc);
% % 
% %     % first, just plot the traces, including the stimulus trace...
% % 
% %     stim = L2_mb(L2_ni_inds(sInd)).cor_stim;
% %     t=[1:length(stim)]/in.xml.framerate; % time in seconds
% % 
% %     act_neurons = active_cells;
% %     Nneurons=length(act_neurons);
% % 
% %     sig=in.dRatio(act_neurons,:);
% %     norm=repmat(mean(sig,2),[1 size(sig,2)]);
% %     sig=sig./norm;
% % 
% %     % Correction
% %     if (length(stim)<size(sig,2))
% %         t = [1:size(sig,2)]/in.xml.framerate;
% %         stim = [stim zeros(1,size(sig,2)-length(stim))];
% %     elseif (length(stim)>size(sig,2))
% %         t = [1:size(sig,2)]/in.xml.framerate;
% %         stim = stim(1:size(sig,2));
% %     end
% % 
% %     figure; hold on;
% %     cm=colormap('lines');
% %     for i=1:Nneurons
% %         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
% %     end
% %     plot(t,stim/3);
% %     ylabel('fractional ratio change');
% %     xlabel('time (s)');
% %     title(sprintf('index =  %d',sInd));
% % 
% %     clear tr;
% %     for i=1:8
% %         tr{i} = aL2ni.rats{sortord(i)}(sInd:(sInd+Nneurons-1),:);
% %     end
% %     t_tr = [0:size(tr{i},2)-1]/100;
% % 
% %     clear tr_n mpt xc yc rc sigang;
% %     for i=1:length(active_cells)
% %         for j=1:8
% %             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
% %         end
% %     end
% % 
% %     mtr=zeros(8,length(t_tr));
% % 
% %     for i=1:length(active_cells)
% %         if (L2_mb(L2_ni_inds(sInd-1+i)).wavelength==490)
% %             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1);
% %             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2);
% %             %         locpos8 = L2_mb(L2_ni_inds(sInd-1+i)).locations;
% %         else
% %             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1)*2;
% %             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2)*2+45;
% %         end
% %     end
% % 
% %     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
% %     x=(x-mean(mean(x))); y=y-mean(mean(y));
% %     x=x/max(max(x))*45*sqrt(2);
% %     y=y/max(max(y))*45*sqrt(2);
% % 
% %     figure; hold on;
% %     plot(xc,yc,'ro-');
% %     for i=1:length(xc)
% %         text(xc(i)+1,yc(i)+1,num2str(active_cells(i)));
% %         circ=2.5*exp(sqrt(-1)*[0:.01:2*pi+0.01]);
% %         plot(xc(i)+real(circ),yc(i)+imag(circ),'k:');
% %     end
% %     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
% %     xlabel('x center');
% %     ylabel('y center');
% %     title(['RF centers: ' num2str(active_cells)]);
% %     set(gca,'dataa',[1 1 1]);
% % 
% %     % get distances
% %     for i = 1:length(active_cells)
% %         cur_center = [xc(i) yc(i)];
% %         all_dists = sqrt((xc - repmat(xc(i),1,length(xc))).^2+...
% %             (yc - repmat(yc(i),1,length(yc))).^2);
% %         all_dists = all_dists(all_dists>0);
% %         dist = min(all_dists);
% %         if(dist>20)
% %             disp(sprintf('error in %d',n));
% %         end
% %         nice_distances = [nice_distances; dist];
% %         disp(sprintf('distance is %0.5g',dist));
% %     end
% % 
% %     figure;
% %     subplot(2,1,1);
% %     imagesc(in.AV1);
% %     axis image;
% %     title('image and masks');
% %     subplot(2,1,2);
% %     maskim=in.AV1;
% %     for i=1:length(active_cells)
% %         maskim(find(in.masks{active_cells(i)}))=(i+1)*((2^12)/(Nneurons+1));
% %     end
% %     hold on;imagesc(flipud(maskim));
% % 
% %     clear mtr;
% % 
% %     %     if(strcmp(input('best?','s'),'y'))
% %     %         best_examples = [best_examples sInd];
% %     %     end
% % 
% %     drawnow;
% %     pause;
% %     close all;
% % end
% 
% vals_mat = [];
% d_vals_mat = [];
% real_diffs = [];
% tt = linspace(-2,2,401);
% for n = 1:length(nice_examples)
% 
%     clc;
%     sInd = nice_examples(n);
%     ii = sInd;
%     fname = aL2ni.name{ii};
%     done = false;
%     active_cells = aL2ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL2ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL2ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL2ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1));
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L2_mb(L2_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
% 
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
% 
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
% %     figure; hold on;
% %     cm=colormap('lines');
% %     for i=1:Nneurons
% %         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
% %     end
% %     plot(t,stim/3);
% %     ylabel('fractional ratio change');
% %     xlabel('time (s)');
% %     title(sprintf('index =  %d',sInd));
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL2ni.rats{sortord(i)}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
% 
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
% 
%     mtr=zeros(8,length(t_tr));
% 
%     order = [1 5 2 6 3 7 4 8];
%     for s = 1:length(active_cells)
%         vals = zeros(8,1);
%         locpos8_temp = L2_mb(L2_ni_inds(sInd-1+s)).locations;
%         locpos8 = zeros(8,1);
%         for i = 1:8
%             [val,locpos8(i)] = min(abs(t_tr-t(locpos8_temp(sortord(i)))));
%         end       
%         for i = 1:8
%             mod_t = t_tr-t_tr(locpos8(i));
%             sig = tr{i}(s,:)-1;
%             zero_ind = find(mod_t==0);
%             ref_ind = find(sig((zero_ind+1):end)>0,1,'first')+zero_ind;
%             vals(i) = mean(sig(ref_ind:min((ref_ind+200),length(sig))))-mean(sig(max((ref_ind-200),1):ref_ind));
% %             range = (ref_ind-200):(ref_ind+200);
% %             figure;plot(mod_t,sig);hold on;
% %             scatter(mod_t(ref_ind),sig(ref_ind));
% %             title(sprintf('angle = %d',round(sortang(i)*180/pi)));
% %             xlim([mod_t(1) mod_t(end)]);
%         end
%         counter = 1;
%         vals_mat = [vals_mat; vals'];
%         d_vals = diff(vals(order));
%         d_vals = d_vals(1:2:end);
%         d_vals_mat = [d_vals_mat; d_vals'];
%         rdiffs = zeros(4,1);
%         for j = 1:2:8
%             sig1 = tr{order(j)}(s,:);
%             mod1_t = t_tr-t_tr(locpos8(order(j)));
%             zero1_ind = find(mod1_t==0);
%             sig2 = tr{order(j+1)}(s,:);
%             mod2_t = t_tr-t_tr(locpos8(order(j+1)));
%             zero2_ind = find(mod2_t==0);
%             ref1_ind = find(sig1((zero1_ind+1):end)>1,1,'first')+zero1_ind;
%             ref2_ind = find(sig2((zero2_ind+1):end)>1,1,'first')+zero2_ind;
%             range1 = max((ref1_ind-200),1):min((ref1_ind+200),length(sig1));
%             range2 = max((ref2_ind-200),1):min((ref2_ind+200),length(sig2));
%             figure;
%             if(length(range1)==401)
%                 plot(tt,sig1(range1),'b');hold on;
%             else
%                 plot(mod1_t(range1),sig1(range1),'b');hold on;
%             end
%             if(length(range2)==401)
%                 plot(tt,sig2(range2),'k');
%             else
%                 plot(mod2_t(range2),sig2(range2),'k');
%             end
%             title([directions_text{counter} ', diff = ',num2str(d_vals(counter))]);
%             legend(angles_text{j},angles_text{j+1});niceaxes;
%             xlabel('time (sec)');ylabel('response (dR/R)');
%             if(strcmp(input('real diff?','s'),'y'))
%                 rdiffs(counter) = 1;
%             end
%             counter = counter+1;
%         end
%         real_diffs = [real_diffs; rdiffs'];
%         close all;
%     end
% end
% 
% save('new_angle_vales','vals_mat','d_vals_mat','real_diffs');
% 
% figure;subplot(2,1,1);bar(vals(order));
% hold on;set(gca,'xTickLabel',angles_text);
% niceaxes;xlim([0.5 8.5]);
% d_vals = diff(vals(order));
% d_vals = d_vals(1:2:end);
% subplot(2,1,2);bar(d_vals);hold on;
% hold on;set(gca,'xTickLabel',directions_text);
% niceaxes;
% 
% %% Plot distribution of distances between RF centers
% 
% % save('angle_vals','vals_mat');
% % save('angle_vals2','vals2_mat');
% % save('angle_vals3','vals3_mat');
% % save('best_traces','best_examples');
% load('L2_distances');
% load('L2_all_distances');
% 
% edges = 0:1:30;
% N = histc(nice_distances,edges);
% figure;bar(edges,N);
% xlim([-5 35]);
% niceaxes;
% title('histogram of angular distance between RF centers: nice examples');
% xlabel('angular distance');xlim([-1 30]);
% ylabel('counts');
% 
% edges = 0:2:30;
% N = histc(all_distances,edges);
% figure;bar(edges,N);
% xlim([-5 35]);
% niceaxes;
% title('histogram of angular distance between RF centers: all');
% xlabel('angular distance');
% ylabel('counts');
% 
% %% Making detailed figures for specific examples
% 
% load('my_examples');
% load('final_MB_curves'); 
% load('best_traces');
% load('angle_vals');
% load('angle_vals2');
% load('angle_vals3');
% 
% % n = 12: nice but RFs are too close to each other
% % n = 15: nice but few (5) cells
% % n = 18: very nice, 6 cells
% % n = 19: very nice, 6 cells (better then prev.)
% % n = 29: 8 cells, not sorted nicely in space
% % n = 34: 10 cells, some are too close together
% % n = 35: 10 cells, separation may not be sufficiently clear from traces
% % n = 37: 10 cells, separation may not be sufficiently clear from traces
% % n = 39: 10 cells, spread isn't perfect, but is good
% % n = 40: 11 cells, spread isn't perfect, but is good
% % n = 41: 11 cells, pretty nice spread
% load('my_examples');
% 
% for n = 19
%     
%     clc;
%     sInd = nice_examples(n);
%     ii = sInd;
%     fname = aL2ni.name{ii};
%     done = false;
%     active_cells = aL2ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL2ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL2ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL2ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1))
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L2_mb(L2_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
% %     figure; hold on;
% %     cm=colormap('lines');
% %     for i=1:Nneurons
% %         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
% %     end
% %     dStim = diff(stim);
% %     inds = find(dStim==1);
% %     text_angs = {'0','180','90','270','315','135','45','225'};
% %     for i = 1:length(inds)
% %         line([t(inds(i)) t(inds(i))],[0 Nneurons+1],'color',[0 0 0],'linewidth',2);
% %         if(i < length(inds))
% %             text(((t(inds(i))+t(inds(i+1)))/2) - 1,Nneurons+0.7,text_angs{i});
% %         else
% %             text((t(inds(i))+t(end)-0.2)/2 - 1,Nneurons+0.7,text_angs{i});
% %         end
% %     end
% %     ylabel('fractional ratio change');
% %     xlabel('time (s)');
% %     ind = find(in.fileloc=='\');
% %     title(in.fileloc((ind(end)+1):end));
% %     ylim([0 Nneurons+1]);
% %     niceaxes;
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL2ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     % Potential supp. figure - how we find the RF center
% %     close all;
% %     for s = 1:length(active_cells)
% %         locpos8_temp = L2_mb(L2_ni_inds(sInd-1+s)).locations;
% %         locpos8 = zeros(8,1);
% %         for i = 1:8
% %             [val,locpos8(i)] = min(abs(t_tr-t(locpos8_temp(i))));
% %         end  
% %         total = locpos8(1)+locpos8(2);
% %         for i=1:2:7
% %             figure; hold on;
% %             plot(t_tr,tr{i}(s,:)','b','lineWidth',2);
% %             plot(t_tr,tr{i+1}(s,end:-1:1)','r','lineWidth',2);
% %             scatter(t_tr(locpos8(i)),tr{i}(s,locpos8(i)),100,[0 0 1],'filled')
% %             scatter(t_tr(total-locpos8(i+1)),tr{i+1}(s,length(t_tr)+locpos8(i+1)-total+1),...
% %                 100,[1 0 0],'filled');
% %             line([t_tr(locpos8(i)) t_tr(locpos8(i))],[0.7 1.3],'color',[0 0 0],'lineWidth',1);
% %             xlabel('+/- time (s)');
% %             ylabel('fractional change');
% %             title(['stim # ' num2str([i i+1]) ', cell = ' num2str(active_cells(s))]);
% %             legend(num2str(ang(i)*180/pi),num2str(ang(i+1)*180/pi));
% %             ylim([0.6 1.4]);
% %         end
% %     end
%     
%     mtr=zeros(8,length(t_tr));
%        
%     for i=1:length(active_cells)
%         if (L2_mb(L2_ni_inds(sInd-1+i)).wavelength==490)
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1);
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2);
% %         locpos8 = L2_mb(L2_ni_inds(sInd-1+i)).locations;
%         else
%             xc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(1)*2;
%             yc(i) = L2_mb(L2_ni_inds(sInd-1+i)).center(2)*2+45;
%         end
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
%     
%     % Create a colored map of ROIs
%     CMask = zeros(in.xml.linesperframe,in.xml.pixperline,3);
%     cm = colormap('lines');
%     nMasks = length(active_cells);
%     for i = 1:nMasks
%         curColor = cm(i,:);
%         curMask = cat(3,curColor(1).*in.masks{active_cells(i)},...
%             curColor(2).*in.masks{active_cells(i)},...
%             curColor(3).*in.masks{active_cells(i)});
%         CMask = CMask + curMask;
%     end
%     figure;
%     imshow(in.AV1,[]);
%     axis image
%     hold on;
%     h = imagesc(CMask./max(CMask(:)));
%     set(h,'AlphaData',0.3);
%     title('image and masks');
%     axis image
% 
%     figure; hold on;
%     plot(xc,yc,'k.-');
%     for i=1:length(xc)
% %         text(xc(i)+3,yc(i),num2str(active_cells(i)));
%         circ = 2.5*exp(sqrt(-1)*[0:10^(-2):(2*pi+10^(-2))]);
%         plot(xc(i)+real(circ),yc(i)+imag(circ),'color',cm(i,:),'linestyle',':','linewidth',2);
%     end
%     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
%     xlabel('x center');
%     ylabel('y center');
% %     title(['RF centers: ' num2str(active_cells)]);
%     title('RF centers');
%     set(gca,'dataa',[1 1 1]);
%     niceaxes;
%     
%     clear mtr;
% 
%     done = false;
% %     close all;
%     while(~done)
%         if(strcmp(input('another direction?','s'),'y'))
%             done = false;
%         else
%             done = true;
%             break;
%         end
%         direction = input('please choose separation direction (horiz,vert,diag1,diag2)','s');
% %         order = input('please choose filter order','s');
% %         cutoff = input('please choose filter cutoff','s');
% %         show_dir_separation3(t_tr,tr,direction,str2num(order),str2num(cutoff));
%         show_dir_separation3(t_tr,tr,direction);
%     end
% 
% %     close all;
% end
% 
% %% Making detailed figures for specific examples - L1_1
% 
% 
% load('my_examples_L1_1'); % nice_examples_L1_1
% nice_examples = nice_examples_L1_1;
% 
% for n = 1:length(nice_examples)
% % for n = 41:41
%     
%     clc;
%     sInd = nice_examples(n);
%     ii = sInd;
%     fname = aL1_1ni.name{ii};
%     done = false;
%     active_cells = aL1_1ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL1_1ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL1_1ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL1_1ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1))
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L1_1_mb(L1_1_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
%     figure; hold on;
%     cm=colormap('lines');
%     for i=1:Nneurons
%         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
%     end
%     dStim = diff(stim);
%     inds = find(dStim==1);
%     text_angs = {'0','180','90','270','315','135','45','225'};
%     for i = 1:length(inds)
%         line([t(inds(i)) t(inds(i))],[0 Nneurons+1],'color',[0 0 0],'linewidth',2);
%         if(i < length(inds))
%             text(((t(inds(i))+t(inds(i+1)))/2) - 1,Nneurons+0.7,text_angs{i});
%         else
%             text((t(inds(i))+t(end)-0.2)/2 - 1,Nneurons+0.7,text_angs{i});
%         end
%     end
%     ylabel('fractional ratio change');
%     xlabel('time (s)');
%     ind = find(in.fileloc=='\');
%     title(in.fileloc((ind(end)+1):end));
%     ylim([0 Nneurons+1]);
%     niceaxes;
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL1_1ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     mtr=zeros(8,length(t_tr));
%        
%     for i=1:length(active_cells)
%         if (L1_1_mb(L1_1_ni_inds(sInd-1+i)).wavelength==490)
%             xc(i) = L1_1_mb(L1_1_ni_inds(sInd-1+i)).centers(1);
%             yc(i) = L1_1_mb(L1_1_ni_inds(sInd-1+i)).centers(2);
% %         locpos8 = L2_mb(L2_ni_inds(sInd-1+i)).locations;
%         else
%             xc(i) = L1_1_mb(L1_1_ni_inds(sInd-1+i)).centers(1)*2;
%             yc(i) = L1_1_mb(L1_1_ni_inds(sInd-1+i)).centers(2)*2+45;
%         end
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
%     
%     % Create a colored map of ROIs
%     CMask = zeros(in.xml.linesperframe,in.xml.pixperline,3);
%     cm = colormap('lines');
%     nMasks = length(active_cells);
%     for i = 1:nMasks
%         curColor = cm(i,:);
%         curMask = cat(3,curColor(1).*in.masks{active_cells(i)},...
%             curColor(2).*in.masks{active_cells(i)},...
%             curColor(3).*in.masks{active_cells(i)});
%         CMask = CMask + curMask;
%     end
%     figure;
%     imshow(in.AV1,[]);
%     axis image
%     hold on;
%     h = imagesc(CMask);
%     set(h,'AlphaData',0.3);
%     title('image and masks');
%     axis image
% 
%     figure; hold on;
%     plot(xc,yc,'k.-');
%     for i=1:length(xc)
% %         text(xc(i)+3,yc(i),num2str(active_cells(i)));
%         circ = 2.5*exp(sqrt(-1)*[0:10^(-2):(2*pi+10^(-2))]);
%         plot(xc(i)+real(circ),yc(i)+imag(circ),'color',cm(i,:),'linestyle',':','linewidth',2);
%     end
%     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
%     xlabel('x center');
%     ylabel('y center');
% %     title(['RF centers: ' num2str(active_cells)]);
%     title('RF centers');
%     set(gca,'dataa',[1 1 1]);
%     niceaxes;
%     
%     clear mtr;
% 
%     done = false;
% %     close all;
%     while(~done)
%         if(strcmp(input('another direction?','s'),'y'))
%             done = false;
%         else
%             done = true;
%             break;
%         end
%         direction = input('please choose separation direction (horiz,vert,diag1,diag2)','s');
% %         order = input('please choose filter order','s');
% %         cutoff = input('please choose filter cutoff','s');
% %         show_dir_separation3(t_tr,tr,direction,str2num(order),str2num(cutoff));
%         show_dir_separation3(t_tr,tr,direction);
%     end
% 
% %     drawnow;
% %     pause;
%     close all;
% end
% 
% %% Making detailed figures for specific examples - L1_5
% 
% 
% load('my_examples_L1_5'); % nice_examples_L1_5
% nice_examples = nice_examples_L1_5;
% 
% % n = 3 - 5 nice cells
% 
% for n = 4:length(nice_examples)
% % for n = 41:41
%     
%     clc;
%     sInd = nice_examples(n);
%     ii = sInd;
%     fname = aL1_5ni.name{ii};
%     done = false;
%     active_cells = aL1_5ni.neuron(ii);
%     while(~done)
%         if(~strcmp(fname,aL1_5ni.name{ii}))
%             done = true;
%         else
%             if(ii < length(aL1_5ni.neuron))
%                 ii = ii+1;
%                 active_cells = [active_cells aL1_5ni.neuron(ii)];
%             else
%                 ii = ii-1;
%                 done = true;
%             end
%         end
%     end
%     active_cells = active_cells(1:(end-1))
% 
%     in=load([pdatapath '/' fname]);
%     in=in.strct;
% 
%     disp(in.fileloc);
% 
%     % first, just plot the traces, including the stimulus trace...
% 
%     stim = L1_5_mb(L1_5_ni_inds(sInd)).cor_stim;
%     t=[1:length(stim)]/in.xml.framerate; % time in seconds
%     
%     act_neurons = active_cells;
%     Nneurons=length(act_neurons);
% 
%     sig=in.dRatio(act_neurons,:);
%     norm=repmat(mean(sig,2),[1 size(sig,2)]);
%     sig=sig./norm;
%     
%     % Correction
%     if (length(stim)<size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = [stim zeros(1,size(sig,2)-length(stim))];
%     elseif (length(stim)>size(sig,2))
%         t = [1:size(sig,2)]/in.xml.framerate;
%         stim = stim(1:size(sig,2));
%     end
% 
%     figure; hold on;
%     cm=colormap('lines');
%     for i=1:Nneurons
%         plot(t,sig(i,:)+i-1,'linewidth',2,'color',cm(i,:));
%     end
%     dStim = diff(stim);
%     inds = find(dStim==1);
%     text_angs = {'0','180','90','270','315','135','45','225'};
%     for i = 1:length(inds)
%         line([t(inds(i)) t(inds(i))],[0 Nneurons+1],'color',[0 0 0],'linewidth',2);
%         if(i < length(inds))
%             text(((t(inds(i))+t(inds(i+1)))/2) - 1,Nneurons+0.7,text_angs{i});
%         else
%             text((t(inds(i))+t(end)-0.2)/2 - 1,Nneurons+0.7,text_angs{i});
%         end
%     end
%     ylabel('fractional ratio change');
%     xlabel('time (s)');
%     ind = find(in.fileloc=='\');
%     title(in.fileloc((ind(end)+1):end));
%     ylim([0 Nneurons+1]);
%     niceaxes;
% 
%     clear tr;
%     for i=1:8
%         tr{i} = aL1_5ni.rats{i}(sInd:(sInd+Nneurons-1),:);
%     end
%     t_tr = [0:size(tr{i},2)-1]/100;
%     
%     clear tr_n mpt xc yc rc sigang;
%     for i=1:length(active_cells)
%         for j=1:8
%             tr_n{i}(j,:)=tr{j}(i,:)-mean(tr{j}(i,:));
%         end
%     end
%     
%     mtr=zeros(8,length(t_tr));
%        
%     for i=1:length(active_cells)
%         if (L1_5_mb(L1_5_ni_inds(sInd-1+i)).wavelength==490)
%             xc(i) = L1_5_mb(L1_5_ni_inds(sInd-1+i)).centers(1);
%             yc(i) = L1_5_mb(L1_5_ni_inds(sInd-1+i)).centers(2);
%         else
%             xc(i) = L1_5_mb(L1_5_ni_inds(sInd-1+i)).centers(1)*2;
%             yc(i) = L1_5_mb(L1_5_ni_inds(sInd-1+i)).centers(2)*2+45;
%         end
%     end
%     
%     [x,y]=meshgrid([1:length(t_tr)],[1:length(t_tr)]);
%     x=(x-mean(mean(x))); y=y-mean(mean(y));
%     x=x/max(max(x))*45*sqrt(2);
%     y=y/max(max(y))*45*sqrt(2);
%     
%     % Create a colored map of ROIs
%     CMask = zeros(in.xml.linesperframe,in.xml.pixperline,3);
%     cm = colormap('lines');
%     nMasks = length(active_cells);
%     for i = 1:nMasks
%         curColor = cm(i,:);
%         curMask = cat(3,curColor(1).*in.masks{active_cells(i)},...
%             curColor(2).*in.masks{active_cells(i)},...
%             curColor(3).*in.masks{active_cells(i)});
%         CMask = CMask + curMask;
%     end
%     figure;
%     imshow(in.AV1,[]);
%     axis image
%     hold on;
%     h = imagesc(CMask);
%     set(h,'AlphaData',0.3);
%     title('image and masks');
%     axis image
% 
%     figure; hold on;
%     plot(xc,yc,'k.-');
%     for i=1:length(xc)
% %         text(xc(i)+3,yc(i),num2str(active_cells(i)));
%         circ = 2.5*exp(sqrt(-1)*[0:10^(-2):(2*pi+10^(-2))]);
%         plot(xc(i)+real(circ),yc(i)+imag(circ),'color',cm(i,:),'linestyle',':','linewidth',2);
%     end
%     plot(45*[-1 -1 1 1 -1],45*[-1 1 1 -1 -1],'k-');
%     xlabel('x center');
%     ylabel('y center');
% %     title(['RF centers: ' num2str(active_cells)]);
%     title('RF centers');
%     set(gca,'dataa',[1 1 1]);
%     niceaxes;
%     
%     clear mtr;
% 
%     done = false;
% %     close all;
%     while(~done)
%         if(strcmp(input('another direction?','s'),'y'))
%             done = false;
%         else
%             done = true;
%             break;
%         end
%         direction = input('please choose separation direction (horiz,vert,diag1,diag2)','s');
% %         order = input('please choose filter order','s');
% %         cutoff = input('please choose filter cutoff','s');
% %         show_dir_separation3(t_tr,tr,direction,str2num(order),str2num(cutoff));
%         show_dir_separation3(t_tr,tr,direction);
%     end
% 
% %     drawnow;
% %     pause;
%     close all;
% end