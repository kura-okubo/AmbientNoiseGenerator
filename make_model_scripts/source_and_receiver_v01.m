%Make ambient noise model for Ambient Noise Generator
%2019 Kurama OKUBO
%%%set environment%%%
clear all;
%clf;
set(0,'DefaultFigureWindowStyle','normal');
%Plot Format
set(0,'DefaultTextFontsize',18, ...
    'DefaultTextFontname','Arial', ...
    'DefaultTextFontWeight','normal', ...
    'DefaultTextFontname','Arial', ...
    'DefaultAxesFontsize',18, ...
    'DefaultAxesFontname','Arial', ...
    'defaultUicontrolFontName','Arial', ...
    'defaultUitableFontName','Arial', ...
    'defaultUipanelFontName','Arial', ...
    'DefaultLineLineWidth', 1.5)

set(0,'defaulttextinterpreter','latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This script makes localized source array with localized receivers, similar
%with Lawrence et al. (2013) Model C in Figure 2. 

%using SI unit (m, s, kg)

%---Model Parameters---%

%Domain size
dim_x = 420e3;
dim_y = 420e3;

%for noise source
NumofSource = 20;
Ls = 50e3; %radius of noise source location (m)
Ls_centre = [150e3, 150e3];

%for receiver
NumofReceiver = 5; %10
Lx = 100e3; %radius of receiver location (m)
Lx_centre = [0, 0];

randomseed_s = 5; %random seed for source locations
randomseed_r = 13; %random seed for receiver locations

problem_name='../EXAMPLE/multi4_source';

SaveFigure = 1;
%----------------------%

fopath = sprintf('./%s/', problem_name);
if isfolder(fopath) == 0; mkdir(fopath); end
if isfolder([fopath, '/inputfiles']) == 0; mkdir([fopath, '/inputfiles']); end

%save const
const.dim_x = dim_x;
const.dim_y = dim_y;
const.NumofSource = NumofSource;
const.Ls = Ls;
const.Ls_centre = Ls_centre;
const.NumofReceiver = NumofReceiver;
const.Lx = Lx;
const.Lx_centre = Lx_centre;

save([fopath,'/inputfiles/const.mat'], 'const');

%calculate location of source and receiver
rng(randomseed_s, 'simdTwister');
randp_s = rand(10000*NumofSource, 2);
count_s = 0;
id_s = 0;

rng(randomseed_r, 'simdTwister');
randp_r = rand(10000*NumofReceiver, 2);
count_r = 0;
id_r = 0;

%minimum distance of sources
mindist_s = 0.15*Ls;
mindist_r = 0.15*Lx;

while count_s < NumofSource
    id_s = id_s + 1; 
    sx_test = Ls * (1 - 2*randp_s(id_s,1))+ Ls_centre(1);
    sy_test = Ls * (1 - 2*randp_s(id_s,2))+ Ls_centre(2);
    if norm([sx_test-+ Ls_centre(1), sy_test-+ Ls_centre(2)], 2) < Ls
        %search the minimum distance        
        flag = 0;
        for i = 1:count_s
            if  norm([sx(i)-sx_test,sy(i)-sy_test], 2) < mindist_s
                flag = 1;
                break;
            end
        end
        
        if flag == 0
            %this location is enough far from any other sources. 
            count_s = count_s + 1;
            sx(count_s) = sx_test;
            sy(count_s) = sy_test;
        end
    end
end

while count_r < NumofReceiver
    id_r = id_r + 1; 
    rx_test = Lx * (1 - 2*randp_r(id_r,1))+ Lx_centre(1);
    ry_test = Lx * (1 - 2*randp_r(id_r,2))+ Lx_centre(2);
    if norm([rx_test-+ Lx_centre(1), ry_test-+ Lx_centre(2)], 2) < Lx
        %search the minimum distance        
        flag = 0;
        for i = 1:count_r
            if  norm([rx(i)-rx_test,ry(i)-ry_test], 2) < mindist_r
                flag = 1;
                break;
            end
        end
        
        if flag == 0
            %this location is enough far from any other sources. 
            count_r = count_r + 1;
            rx(count_r) = rx_test;
            ry(count_r) = ry_test;
        end
    end
end
%%

%analyze the distribution of distance between sources
scount = 1;
sdist=[]
for i = 1:NumofSource
    for j = i:NumofSource
        if i==j 
            continue;
        else
            sdist(scount) = norm([sx(i)-sx(j),sy(i)-sy(j)], 2);
            scount = scount + 1;
        end
    end
end

rcount = 1;
for i = 1:NumofReceiver
    for j = i:NumofReceiver
        if i==j 
            continue;
        else
            rdist(rcount) = norm([rx(i)-rx(j),ry(i)-ry(j)], 2);
            rcount = rcount + 1;
        end
    end
end

%%
%compute array response function (ARF) for vertical signal ([kx, ky] = [0, 0])
gridnum = 101;
f_ARF = 0.1;
maxslowness = 0.1*1e-3; %[s/m]
kmax = 2*pi*f_ARF*maxslowness;
kx = linspace(-kmax, kmax, gridnum); 
ky = linspace(-kmax, kmax, gridnum); 
E = zeros(gridnum,gridnum);

%Convert to dB
A0 = 1; %because |k0| = 0

icount = 1;
for i = 1:gridnum    
    for j = 1:gridnum
        Etemp = 0;
        for m = 1:NumofReceiver
            kv = [kx(i), ky(j)];
            rv = [rx(m)-Lx_centre(1), ry(m)-Lx_centre(2)];
            Etemp = Etemp + exp(2*pi*1i*dot(kv, rv));
        end
        
        Etemp = Etemp*conj(Etemp)/NumofReceiver^2;
        E(i,j) = 10*log10(Etemp/A0);
        icount = icount+1;
    end
end


fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.5 0.5];
clf(fig,'reset'); cla(fig,'reset'); hold on;
axis square;


[pkx, pky] = meshgrid(kx, ky);

contourf(pkx*1e3, pky*1e3, E);
ax1 = gca;


XLimit = [-kmax, kmax];
YLimit = [-kmax, kmax];

ax1.XLim = XLimit*1e3;
ax1.YLim = YLimit*1e3;

xlabel('kx (1/km)');
ylabel('ky (1/km)');
title(sprintf("Array Response Function for %4.2f [Hz]", f_ARF));
box on;

%superpose polar coordinate
hold on;
ac = [100 100 100]/255;

t = 0 : .01 : 2 * pi;
hp = polar(t, (kmax) * ones(size(t)) *1e3, 'k');
hp = polar(t, 0.33333*(kmax) * ones(size(t)) *1e3, 'k');
hp = polar(t, 0.6666*(kmax) * ones(size(t)) *1e3, 'k');
text(0, 0.33333*(kmax)*1e3, sprintf('%4.4f [s/km]', 0.33333*(kmax)*1e3), 'Color', 'g');
text(0, 0.66666*(kmax)*1e3, sprintf('%4.4f [s/km]', 0.66666*(kmax)*1e3), 'Color', 'g');
text(0, 1.0*(kmax)*1e3, sprintf('%4.4f [s/km]', 1.0*(kmax)*1e3), 'Color', 'g');

ptheta = linspace(0, pi, 7);
pr = kmax*ones(length(ptheta), 1)';

for i = 1:length(ptheta)
    rr = pr(i);
    tt = ptheta(i);
    plot(rr*[-cos(tt), cos(tt)]*1e3, rr*[-sin(tt), sin(tt)]*1e3, 'Color', ac,...
        'LineWidth', 1.0);
end

c1 = colorbar;
colormap(jet);
caxis([-24, 0]);
ylabel(c1, 'Power [dB]')

if (SaveFigure)
    figdir_depth = sprintf('%s/model_fig/png', fopath);
    fodir = [figdir_depth,'/'];
    if isfolder(fodir) == 0; mkdir(fodir); end
    set(gcf, 'Color', 'w');
    foname = sprintf('%s/ARF.%s', fodir, 'png');
    export_fig(foname,'-r200');
end

%%
fig = figure(2);
fig.Units = 'normalized';
fig.Position = [0 1 0.6 0.6];
clf(fig,'reset'); cla(fig,'reset'); hold on;
axis square

h = scatter(sx/1e3, sy/1e3, 3, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'w');
hold on;
h = scatter(rx/1e3, ry/1e3, 30, 'v', 'MarkerEdgeColor', 'b');

%plot circle line
xc = Ls * cos(linspace(0,2*pi,100))+ Ls_centre(1);
yc = Ls * sin(linspace(0,2*pi,100))+ Ls_centre(2);
h = plot(xc/1e3, yc/1e3, 'k-', 'LineWidth',1.0);

xc = Lx * cos(linspace(0,2*pi,100))+ Lx_centre(1);
yc = Lx * sin(linspace(0,2*pi,100))+ Lx_centre(2);
h = plot(xc/1e3, yc/1e3, 'k-', 'LineWidth',1.0);

ax1 = gca;
XLimit = [-0.5*dim_x, 0.5*dim_x]/1e3;
YLimit = [-0.5*dim_y, 0.5*dim_y]/1e3;

ax1.XLim = XLimit;
ax1.YLim = YLimit;

xlabel('x (km)');
ylabel('y (km)');

box on;

if (SaveFigure)
    figdir_depth = sprintf('%s/model_fig/png', fopath);
    fodir = [figdir_depth,'/'];
    if isfolder(fodir) == 0; mkdir(fodir); end
    set(gcf, 'Color', 'w');
    foname = sprintf('%s/Source_and_Receiver_loc.%s', fodir, 'png');
    export_fig(foname,'-r200');
end

fig = figure(3);
fig.Units = 'normalized';
fig.Position = [0 1 0.4 0.4];
clf(fig,'reset'); cla(fig,'reset'); hold on;
axis square
subplot(2,1,1)
p = histogram(sdist/1e3, 'Normalization','count');
xlabel('distance between sources (km)');
ylabel('count');

subplot(2,1,2)
p = histogram(rdist/1e3, 'Normalization','count');
xlabel('distance between receivers (km)');
ylabel('count');

if (SaveFigure)
    figdir_depth = sprintf('%s/model_fig/png', fopath);
    fodir = [figdir_depth,'/'];
    if isfolder(fodir) == 0; mkdir(fodir); end
    set(gcf, 'Color', 'w');
    foname = sprintf('%s/Source_receiver_dist_distribution.%s', fodir, 'png');
    export_fig(foname,'-r200');
end


%%
%output source locations
fo_source = [fopath+"/inputfiles/source_loc.in"];
if isfolder(fopath) == 0; mkdir(fopath); end

fileID = fopen(fo_source,'w');
fprintf(fileID,'#x[m] y[m]\n');
for i = 1:length(sx)
    fprintf(fileID, sprintf('%12.8f, %12.8f\n', sx(i), sy(i)));
end
fclose(fileID);

%output receiver locations
fo_receiver= [fopath+"/inputfiles/receiver_loc.in"];
if isfolder(fopath) == 0; mkdir(fopath); end

fileID = fopen(fo_receiver, 'w');
fprintf(fileID,'#x[m] y[m]\n');
for i = 1:length(rx)
    fprintf(fileID, sprintf('%12.8f, %12.8f\n', rx(i), ry(i)));
end
fclose(fileID);

