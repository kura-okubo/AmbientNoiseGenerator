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

problem_name = "multi3_source";

Time_ID = [1];
Day_ID = [1:2];
Unit_ID = [1:5];

%for plotting
SaveFigure = 1;

f_inputdir = "../EXAMPLE/"+problem_name+"/OUTPUT_FILES/";

Figdir = "../EXAMPLE/"+problem_name+"/OUTPUT_FILES/figs";
FileFormat = 'png';

fiR = f_inputdir+"Receiver.h5";
fiCC1 = f_inputdir+"CC1.h5";
fiCC2 = f_inputdir+"CC2.h5";
fiCC3 = f_inputdir+"CC3.h5";

%Load data

t =  hdf5read(fiR, "t");

ReceiverID = [1:5];
TimeID = 1;
DayID = 1;
UnitID = 3;

u = zeros(length(t), length(ReceiverID));

for i = ReceiverID
    
    dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/U.%02d/U", TimeID, DayID, UnitID, i);
    u(:, i) =  hdf5read(fiR, dataname);
end


t_cc =  hdf5read(fiCC3, "Lag_time");
CC3dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/CC3.V%02d-%02d-%02d/CC1_v1_all", TimeID, DayID, UnitID, 5, 3, 4);
cc_v1_all = hdf5read(fiCC3, CC3dataname);

CC3dataname = sprintf("TimeID%02d/Day%04d/UnitID%04d/CC3.V%02d-%02d-%02d/CC1_v1_windowed", TimeID, DayID, UnitID, 5, 3, 4);
cc_v1_windowed = hdf5read(fiCC3, CC3dataname);

C1_Normalize = 8e-5;

fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.5 0.3];
clf(fig,'reset'); cla(fig,'reset'); hold on;

plot(t_cc/60, cc_v1_all./C1_Normalize, "k-");
plot(t_cc/60, cc_v1_windowed./C1_Normalize, "r-");


box on;
ax1 = gca;
XLimit = [-1200, 1200]/60;
YLimit = [-1, 1];

ax1.XLim = XLimit;
ax1.YLim = YLimit;

xlabel('Lag time (min)');
ylabel('Coherency');


%%
%plot noise data

fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.5 0.6];
clf(fig,'reset'); cla(fig,'reset'); hold on;

u_normalized = max(max(u));
trace_span = 1;
for i = ReceiverID
    plot(t/60, 0.8 * u(:,i)./u_normalized + (i) * trace_span, 'k-');
end

box on;
ax1 = gca;
XLimit = [0, 60];
YLimit = [-1, 1];

ax1.XLim = XLimit;
%ax1.YLim = YLimit;

xlabel('Time (min)');
ylabel('Receiver ID');

if (SaveFigure)
    fodir = [Figdir+'/'];
    mkdir(fodir);
    set(gcf, 'Color', 'w');
    foname = sprintf('%s/noisefield_TimeID%02d_Day%04d_UnitID%04d.%s',fodir,...
    TimeID, DayID, UnitID, 'png');
    export_fig(foname,'-r200');
end

%%

fig = figure(2);
fig.Units = 'normalized';
fig.Position = [0 1 0.5 0.3];
clf(fig,'reset'); cla(fig,'reset'); hold on;

%plot spectrum
downsamplerate = 2;
dt = downsamplerate * (t(2) - t(1));
fs = 1/dt;
U_ID=2;

spectrogram(u(1:downsamplerate:end, U_ID)./max(abs(u(:, U_ID))) ,kaiser(60,16), 20, 128, fs, 'yaxis');

ax = gca;

XLimit = [0, 60];
YLimit = [0.1, 2];
ax.YScale = 'log';
ax.XLim = XLimit;
ax.YLim = YLimit;
xlabel('Time (min)');
ylabel('Frequency (Hz)');
caxis([-40, 0]);

yt = get(gca, 'YTick');
%set(gca, 'YTick',yt, 'YTickLabel',yt*1E-3);

if (SaveFigure)
    fodir = [Figdir+'/'];
    mkdir(fodir);
    set(gcf, 'Color', 'w');
    foname = sprintf('%s/noisespectrogram_TimeID%02d_Day%04d_UnitID%04d.%s',fodir,...
    TimeID, DayID, U_ID, 'png');
    export_fig(foname,'-r200');
end
