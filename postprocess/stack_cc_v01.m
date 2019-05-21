%Stack cross-correlation function
%2019.05.20 Kurama OKUBO

%%%set environment%%%
clear all;
%clf;
set(0,'DefaultFigureWindowStyle','normal');
set(0,'defaulttextinterpreter','latex');

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

set(0,'defaulttextinterpreter','tex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

problem_name = "multi4_source";

Time_ID = [1];
Day_ID = [1:2];
Unit_ID = [1:12];

%for plotting
SaveFigure = 1;


f_inputdir = "../EXAMPLE/"+problem_name+"/OUTPUT_FILES/";

Figdir = "../EXAMPLE/"+problem_name+"/OUTPUT_FILES/figs";
FileFormat = 'png';

fiCC1 = f_inputdir+"CC1.h5";
fiCC2 = f_inputdir+"CC2.h5";
fiCC3 = f_inputdir+"CC3.h5";

%Load data
CC1_lagtime = hdf5read(fiCC1, "Lag_time");
CC1_Receiver_pair = hdf5read(fiCC1, "Receiver_pair");

CC3_lagtime = hdf5read(fiCC3, "Lag_time");
CC3_Receiver_group = hdf5read(fiCC3, "Receiver_group");

%Target receiver pair
R1 = 3;
R2 = 4;

%Target virtual source
V1 = [5];

%search ID to read hdf5 dataset
for i = 1:length(CC1_Receiver_pair)
    if (CC1_Receiver_pair(i,1)==R1 && CC1_Receiver_pair(i,2)==R2) || (CC1_Receiver_pair(i,1)==R2 && CC1_Receiver_pair(i,2)==R1)
        PairID = i;
    end
end

GroupIDcount = 0;
for i = 1:length(CC3_Receiver_group)
    if ismember(CC3_Receiver_group(i,1), V1)
        if (CC3_Receiver_group(i,2)==R1 && CC3_Receiver_group(i,3)==R2) || (CC3_Receiver_group(i,2)==R2 && CC3_Receiver_group(i,3)==R1)
            GroupIDcount = GroupIDcount + 1;
            GroupID(GroupIDcount) = i;
        end
    end
end

%%
%stack cc1 and cc3
NumofStack = max(Time_ID) * max(Day_ID) * max(Unit_ID);
CC1stack = zeros(length(CC1_lagtime), NumofStack);
CC3stack = zeros(length(CC3_lagtime), NumofStack);

stackcount = 0;
for i = Time_ID
    for j = Day_ID
        for k = Unit_ID
            
            %For C1
            if CC1_Receiver_pair(PairID, 1) == R1
                CC1dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/CC1.%02d-%02d/CC1", i, j, k, R1, R2);
                CC1_temp = hdf5read(fiCC1, CC1dataname);
                CC1attrloc = sprintf("/TimeID%02d/Day%04d/UnitID%04d/CC1.%02d-%02d", i, j, k, R1, R2);
                r1x =   h5readatt(fiCC1, CC1attrloc, "r1x");
                r1y =   h5readatt(fiCC1, CC1attrloc, "r1y");
                r2x =   h5readatt(fiCC1, CC1attrloc, "r2x");
                r2y =   h5readatt(fiCC1, CC1attrloc, "r2y");
                dist = norm([(r1x-r2x), (r1y-r2y)], 2);
            else
                CC1dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/CC1.%02d-%02d/CC1", i, j, k, R2, R1);
                CC1_temp = fliplr(hdf5read(fiCC1, CC1dataname));
                CC1attrloc = sprintf("/TimeID%02d/Day%04d/UnitID%04d/CC1.%02d-%02d", i, j, k, R2, R1);
                r1x =   h5readatt(fiCC1, CC1attrloc, "r1x");
                r1y =   h5readatt(fiCC1, CC1attrloc, "r1y");
                r2x =   h5readatt(fiCC1, CC1attrloc, "r2x");
                r2y =   h5readatt(fiCC1, CC1attrloc, "r2y");
                dist = norm([(r1x-r2x), (r1y-r2y)], 2);
            end
            

            %For C3
            for l = 1:GroupIDcount
                if CC3_Receiver_group(GroupID(l), 2) == R1
                    CC3dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/CC3.V%02d-%02d-%02d/CC3", i, j, k, V1, R1, R2);
                    CC3_temp(:, l) = hdf5read(fiCC3, CC3dataname);
                else
                    CC3dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/CC3.V%02d-%02d-%02d/CC3", i, j, k, V1, R2, R1);
                    CC3_temp(:, l) = fliplr(hdf5read(fiCC3, CC3dataname));
                end
            end
            
            %Stack C1 and C3
            stackcount = stackcount + 1;
            CC1stack(:, stackcount) = CC1_temp;
            
            for l = 1:GroupIDcount
                CC3stack(:, stackcount) = CC3stack(:, stackcount) + CC3_temp(:, l);
            end
            
        end
    end
end

%%
%Plot C1 and C3
            
fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.5 0.6];
clf(fig,'reset'); cla(fig,'reset'); hold on;
axis square


%plot C1

%cmap = jet(NumofStack);

colorname = "pararainbow.json";
colorname = "pararainbow_zerowhite.json";
cmap = importColormapFromParaview(colorname, 0, 0, 0, NumofStack);
cmap(1,:) = [0.8, 0.8, 0.8];

%C1_Normalize = 1.0e-04;
%C3_Normalize = 6e-08;

C1_Normalize = 0.05;
C3_Normalize = 0.075;

for i = 1:NumofStack
    
    StackNum = double(i);
    
    subplot(2,1,1)
    hold on;
    plot(CC1_lagtime, sum(CC1stack(:, 1:i), 2)./ C1_Normalize ./StackNum, '-', 'Color', cmap(i, :));

    ax1 = gca;
    XLimit = [-60, 60];
    YLimit = [-1, 1];

    ax1.XLim = XLimit;
    ax1.YLim = YLimit;

    xlabel('Lag time (s)');
    ylabel('Normalized Coherency');
    title(sprintf("C1: R %02d-%02d Dist:%4.2f [km] Stack: %02d", R1, R2, dist/1e3, i))
    box on;

    %plot C3
    subplot(2,1,2)
    hold on;
    plot(CC3_lagtime, sum(CC3stack(:, 1:i), 2)./C3_Normalize./StackNum, '-', 'Color', cmap(i, :));

    ax1 = gca;
    ax1.XLim = XLimit;
    ax1.YLim = YLimit;

    xlabel('Lag time (s)');
    ylabel('Normalized Coherency');
    title(sprintf("C3: V %02d R %02d-%02d Dist:%4.2f [km] Stack: %02d", V1, R1, R2, dist/1e3, i))
    box on;
    
    h = colorbar('eastoutside');
    colormap(importColormapFromParaview(colorname, 0, 0, 0, 201));
    clear caxis;
    caxis([1 NumofStack]);
    hp = h.Position;
    h.Position = [hp(1)+0.05, hp(2) + 0.15, hp(3), hp(4)+0.2];
    ylabel(h, 'Stack Number (hourly stack)');
    set(h,'YTick',[0:1:NumofStack]);
    
    if SaveFigure
        fodir = [Figdir+'/'];
        if isfolder(fodir) == 0; mkdir(fodir); end
        set(gcf, 'Color', 'w');
        foname = sprintf('%s/C1andC3_stack_%04d.%s', fodir, i, 'png');
        export_fig(foname,'-r200');
    end
    pause(0.1)

end

            
%%
datasetpath = '../dataset/finitefaultenergy/';
rupturetype = './';
casetype = './';
depthIDList = [2]; %km



if (SaveFigure)
    figdir_depth = sprintf('../fig/STF/png');
    fodir = [figdir_depth,'/'];
    mkdir(fodir);
    set(gcf, 'Color', 'w');
    %print('-painters', '-depsc', [fodir, sprintf('/EnergyBalance_barchart_Timestep%02d_continuum.%s',...
    %    TimeID, 'eps')]);
    foname = sprintf('%s/STF_finitefault.%s', fodir, 'eps');
    export_fig(foname,'-r200');
end
