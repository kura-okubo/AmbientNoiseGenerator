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

%This script makes time-dependent

%using SI unit (m, s, kg)

%---Model Parameters---%

Time_ID = [1, 2]; %ID of time period for changing v over time;
init_day = [1, 5]; %init julian day of each time period
end_day = [5, 10]; %end julian day of each time period

%velocity%
dv_V=-0.02; %for dv/V test; true value

%v0 = [3333, (1+dv_V)*3333]; %[m/s]
v0 = [3333, 3333]; %[m/s]
IsDispersion = 1;
Dispersion_tickness = [10e3, 10e3]; %[m]

%scatter
scatter_pattern = "random";
scatter_randomseed = [12, 12];
scatter_num = 60;
scatter_mindist= 8e3; % [m]minimum distance between scattering point
scatter_meanstrength = [4e7, 4e7]; %[6e7, 6e7];
scatter_variance = 0.1 .* scatter_meanstrength;

%attenuation
IsAttenuation = 1;
Q = [500, 250]; %[250, 250];

problem_name='../EXAMPLE/coda_test';

SaveFigure = 1;

%----------------------%

%load input
fopath = sprintf('./%s/', problem_name);
if isfolder(fopath) == 0; error("run source_and_receiver.m first to make problem."); end
if isfolder([fopath, '/inputfiles']) == 0; error("run source_and_receiver.m first to make problem."); end

load([fopath, '/inputfiles/const.mat']);
A = importdata([fopath, '/inputfiles/source_loc.in'], ',', 1);
sx = A.data(:, 1);
sy = A.data(:, 2);

A = importdata([fopath, '/inputfiles/receiver_loc.in'], ',', 1);
rx = A.data(:, 1);
ry = A.data(:, 2);

%loop timeid
for l = 1:length(Time_ID)
%for l = 1
    timeid = Time_ID(l);
    
    fofoldername = sprintf('%s/inputfiles/TimeID_%02d', fopath, timeid);
    
    if isfolder(fofoldername) == 0; mkdir(fofoldername); end
    
    %velocity structure
    D = importdata('./dispersioncurve/dispersion_profile.in', ',', 1);
    
    if IsDispersion
        period = D.data(:,1) .* Dispersion_tickness(l) ./v0(l);
        vel = D.data(:,2) .* v0(l);
    else
        period = D.data(:,1);
        vel = ones(length(period),1) .* 0.919402*v0(l);
    end
    
    %scatter distribution
    
    %define scatter distribution area
    %assume ellipse
    t=0:10:360;
    n=length(t);
    
    %arbitrary chose size of ellipse
    %a=1.2*norm([const.Ls_centre(1)-const.Lx_centre(1), const.Ls_centre(2)-const.Lx_centre(2)], 2);
    a=2.0*norm([const.Ls_centre(1)-const.Lx_centre(1), const.Ls_centre(2)-const.Lx_centre(2)], 2);
    b=0.6*a;
    theta=atand((const.Ls_centre(1)-const.Lx_centre(1))/(const.Ls_centre(2)-const.Lx_centre(2)));
    
    x0 = 0.8*(const.Ls_centre(1)+const.Lx_centre(1))/2;
    y0 = 0.8*(const.Ls_centre(2)+const.Lx_centre(2))/2;
    c=[x0*ones(1,n);y0*ones(1,n)];
    x1=[a*cosd(t);
        b*sind(t)];
    
    R=[cosd(theta) -sind(theta);
        sind(theta) cosd(theta)];
    X=R*x1+c;
    
    rng(scatter_randomseed(l), 'simdTwister');
    randp_s = rand(100*scatter_num, 3);
    sc_strength = normrnd(scatter_meanstrength(l),scatter_variance(l), scatter_num, 1);
    sc_strength(sc_strength < 0) = 0.0;
    count_s = 0;
    id_s = 0;
    
    while count_s < scatter_num
        id_s = id_s + 1;
        scx_test = a * (1 - 2*randp_s(id_s,1))+ x0;
        scy_test = a * (1 - 2*randp_s(id_s,2))+ y0;
        if inpolygon(scx_test, scy_test, X(1,:), X(2,:))
            %search the minimum distance
            flag = 0;
            for i = 1:count_s
                if  norm([scx(i)-scx_test,scy(i)-scy_test], 2) < scatter_mindist
                    flag = 1;
                    break;
                end
            end
            
            if flag == 0
                %this location is enough far from any other sources.
                count_s = count_s + 1;
                scx(count_s) = scx_test;
                scy(count_s) = scy_test;
            end
        end
    end
    
    %attenuation
    velintp = @(x) interp1(period, vel, x, 'spline');
    if IsAttenuation
        Q_period = linspace(min(period), max(period), 2*length(period)) ;
        omega = 2*pi./Q_period;
        Q_alpha = omega ./ (2*velintp(omega)*Q(l)); %[1/m]
    else
        Q_period = linspace(min(period), max(period), 2*length(period)) ;
        omega = 2*pi./Q_period;
        %put huge Q to be none
        Q_alpha = omega ./ (2*velintp(omega)*1e6); %[1/m]
    end
    
    %save input files
    fo_model_sc = sprintf('%s/velocity.in', fofoldername);
    fileID = fopen(fo_model_sc, 'w');
    fprintf(fileID,'#period[s] phase velocity[m/s]\n');
    for i = 1:length(period)
        fprintf(fileID, sprintf('%12.8f, %12.8f\n', period(i), vel(i)));
    end
    fclose(fileID);
    
    fo_model_sc = sprintf('%s/scatter_loc.in', fofoldername);
    fileID = fopen(fo_model_sc, 'w');
    fprintf(fileID,'#x[m] y[m] scatter_strength\n');
    for i = 1:length(scx)
        fprintf(fileID, sprintf('%12.8f, %12.8f, %12.8f\n', scx(i), scy(i), sc_strength(i)));
    end
    fclose(fileID);
    
    fo_model_sc = sprintf('%s/attenuation.in', fofoldername);
    fileID = fopen(fo_model_sc, 'w');
    fprintf(fileID,'#period[s] alpha[1/m]\n');
    for i = 1:length(Q_period)
        fprintf(fileID, sprintf('%12.8e, %12.8e\n', Q_period(i), Q_alpha(i)));
    end
    fclose(fileID);
    
    %plot figure
    if SaveFigure
        %plot model
        fig = figure(1);
        fig.Units = 'normalized';
        fig.Position = [0 1 0.6 0.6];
        clf(fig,'reset'); cla(fig,'reset'); hold on;
        axis square
        
        %plot scatter points
        scatter(scx/1e3, scy/1e3, 30, sc_strength,'filled',...
            'LineWidth',0.5, 'MarkerEdgeColor', 'none');
        cmap = colormap('gray');
        colormap(flipud(cmap));
        
        h=colorbar;
        ylabel(h, 'Scatter strength');
        
        scatter(sx/1e3, sy/1e3, 60, 'p',...
            'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
        
        scatter(rx/1e3, ry/1e3, 100, 'v',...
            'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
        
        %plot circle line
        xc = const.Ls * cos(linspace(0,2*pi,100))+ const.Ls_centre(1);
        yc = const.Ls * sin(linspace(0,2*pi,100))+ const.Ls_centre(2);
        plot(xc/1e3, yc/1e3, 'k-', 'LineWidth',1.0);
        
        xc = const.Lx * cos(linspace(0,2*pi,100))+ const.Lx_centre(1);
        yc = const.Lx * sin(linspace(0,2*pi,100))+ const.Lx_centre(2);
        plot(xc/1e3, yc/1e3, 'k-', 'LineWidth',1.0);
        
        plot(X(1,:)/1e3,X(2,:)/1e3, 'k-', 'LineWidth',1.0);
        
        ax1 = gca;
        XLimit = 1.5*[-0.5*const.dim_x, 0.5*const.dim_x]/1e3;
        YLimit = 1.5*[-0.5*const.dim_y, 0.5*const.dim_y]/1e3;
        
        ax1.XLim = XLimit;
        ax1.YLim = YLimit;
        
        xlabel('x (km)');
        ylabel('y (km)');
        
        box on;
        
        figdir_depth = sprintf('%s/model_fig/', fofoldername);
        fodir = [figdir_depth,'/'];
        if isfolder(fodir) == 0; mkdir(fodir); end
        set(gcf, 'Color', 'w');
        foname = sprintf('%s/structuremodel.%s', fodir, 'png');
        export_fig(foname,'-r200');
        
        fig = figure(2);
        fig.Units = 'normalized';
        fig.Position = [0 1 0.6 0.6];
        clf(fig,'reset'); cla(fig,'reset'); hold on;
        subplot(2,1,1)
        plot(period, vel/1e3, 'k-');
        xlabel('period (s)');
        ylabel('velocity (km/s)');
        
        subplot(2,1,2)
        plot(Q_period, Q_alpha*1e3, 'r-');
        xlabel('period (s)');
        ylabel('$\alpha$ (1/km)');
        
        set(gcf, 'Color', 'w');
        foname = sprintf('%s/velosityandattenuation_dispersioncurve.%s', fodir, 'png');
        export_fig(foname,'-r200');

        
    end
    
end

