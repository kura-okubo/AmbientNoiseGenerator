%Make ambient noise config for Ambient Noise Generator
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

%This script makes configure file for ambient noise modeling

%---Model Parameters---%
problem_name='../EXAMPLE/multi4_source';

%input_filenames
in_source = "source_loc.in";
in_receiver = "receiver_loc.in";
in_velocity = "velocity.in";
in_scatter = "scatter_loc.in";
in_attenuation = "attenuation.in";

%Time index infomation
Time_ID = [1, 2]; %ID of time period for changing v over time;
init_day = [1, 6]; %init julian day of each time period
end_day = [5, 10]; %end julian day of each time period

%density to be here
rho = [3000, 3000];

%sampling information
sampling_frequency = 10; %[Hz]
minimum_frequency = 0.02; %[Hz]

%attenuation option
IsAttenuation = "true";

%Normalization by spectral amplitude |u1||u2|
IsSpectralNormalization = "true";

%source information
Sourcetype = "ricker";

IsGaussianNoise = "false"; %apply gaussian noise

Source_randomseed = [10, 11];

Source_amp_mean = 10;
Source_amp_variance = 0.01*Source_amp_mean;

Source_peakfreq_mean = 0.5; %[Hz]
Source_peakfreq_variance = 20; %[percent] of mean peak freq 

Stacking_Frequency = 1; % unit of stacking [hours]
Source_average_num_per_unithour = 30; %source emit this num per Stacking_Frequency [1/unitstackfreq].

%C1 and C3 option
Synthesize_unit_duration = 1*60*60; %[s] the noise is synthesized every this unit time.
T_for_gf = 300; %[s] maximum duration for a green's function between source and receiver. 
                %This gf is merged over Synthesize_unit_duration with
                %random activation time

C3_truncate_alpha = 2.0; %trancate ballistic wave with alpha * r_distance / mean(vel)
C3_max_time_lag = 15*60; %[s] max time lag for C3

%Output option
IsSaveRawTrace  = "true"; %save raw traces
IsSaveCC        = "true"; %save cross_correlation_traces

IsSaveFigure = 0; %for checking this configration
%----------------------%

fopath = sprintf('./%s/', problem_name);
if isfolder(fopath) == 0; error("run source_and_receiver.m first to make problem."); end
if isfolder([fopath, '/inputfiles']) == 0; error("run source_and_receiver.m first to make problem."); end

load([fopath, '/inputfiles/const.mat']);

fodir = sprintf("%s/model_fig/png", problem_name);

%check Time continuous
for i = 1:length(Time_ID)-1
    if init_day(i+1) ~= end_day(i)+1
        error('init_day and end_day should be continuous');
    end
end

%check source config
rng(Source_randomseed(1), 'simdTwister');

fig = figure(1);
fig.Units = 'normalized';
fig.Position = [0 1 0.6 0.6];
clf(fig,'reset'); cla(fig,'reset'); hold on;

subplot(2,1,1)
histogram(normrnd(Source_amp_mean, Source_amp_variance, const.NumofSource));
xlabel('Source amplitude');
ylabel('count');

subplot(2,1,2)

mu = log((Source_peakfreq_mean^2)/sqrt(Source_peakfreq_variance+Source_peakfreq_mean^2));
sigma = sqrt(log(Source_peakfreq_variance/(Source_peakfreq_mean^2)+1));

sourcefreqdist = lognrnd(mu, sigma, const.NumofSource);

[~,edges] = histcounts(log(sourcefreqdist));
histogram(sourcefreqdist, exp(edges))
 
ax = gca;
ax.XScale = 'log';
ax.XLim = [0.01, 5.0];
xlabel('Source peak frequency');
ylabel('count');

if IsSaveFigure
    figdir_depth = sprintf('%s', fodir);
    fodir = [figdir_depth,'/'];
    if isfolder(fodir) == 0; mkdir(fodir); end
    set(gcf, 'Color', 'w');
    foname = sprintf('%s/source_statistics.%s', fodir, 'png');
    export_fig(foname,'-r200');
end

%Save config input file
fo_config = [fopath, '/inputfiles/config.in'];
fileID = fopen(fo_config, 'w');
fprintf(fileID,'#config parameters made at %s\n', string(datetime));
fprintf(fileID,'#Input file names\n');
fprintf(fileID, sprintf('in_source              =%s\n',in_source));
fprintf(fileID, sprintf('in_receiver            =%s\n',in_receiver));
fprintf(fileID, sprintf('in_velocity            =%s\n',in_velocity));
fprintf(fileID, sprintf('in_scatter             =%s\n',in_scatter));
fprintf(fileID, sprintf('in_attenuation         =%s\n',in_attenuation));
fprintf(fileID,'#Simulation time informations\n');
fprintf(fileID, sprintf('NumofTime_ID           =%s\n',(num2str(length(Time_ID)))));
fprintf(fileID, sprintf('Time_ID                =%s\n',regexprep(num2str(Time_ID),'\s+',' ')));
fprintf(fileID, sprintf('init_day               =%s\n',regexprep(num2str(init_day),'\s+',' ')));
fprintf(fileID, sprintf('end_day                =%s\n',regexprep(num2str(end_day),'\s+',' ')));
fprintf(fileID, sprintf('rho     =%s\n',regexprep(num2str(rho),'\s+',' ')));
fprintf(fileID,'#Sampling frequency\n');
fprintf(fileID, sprintf('sampling_frequency     =%8.4f\n',sampling_frequency));
fprintf(fileID, sprintf('minimnum_frequency     =%8.4f\n',minimum_frequency));
fprintf(fileID, sprintf('Synthesize_unit_duration =%8.4f\n', Synthesize_unit_duration));
fprintf(fileID, sprintf('T_for_gf =%8.4f\n', T_for_gf));
fprintf(fileID,'#Attenuation option\n');
fprintf(fileID, sprintf('IsAttenuation          =%s\n',IsAttenuation));
fprintf(fileID, sprintf('IsSpectralNormalization          =%s\n',IsSpectralNormalization));
fprintf(fileID,'#Gaussian noise option\n');
fprintf(fileID, sprintf('IsGaussianNoise        =%s\n',IsGaussianNoise));
fprintf(fileID,'#Source option\n');
fprintf(fileID, sprintf('Sourcetype             =%s\n',Sourcetype));
fprintf(fileID, sprintf('Source_randomseed      =%s\n', regexprep(num2str(Source_randomseed),'\s+',' ')));
fprintf(fileID, sprintf('Source_amp_mean        =%8.4f\n',Source_amp_mean));
fprintf(fileID, sprintf('Source_amp_variance    =%8.4f\n',Source_amp_variance));
fprintf(fileID, sprintf('Source_peakfreq_mean   =%8.4f\n',Source_peakfreq_mean));
fprintf(fileID, sprintf('Source_peakfreq_variance=%8.4f\n',Source_peakfreq_variance));
fprintf(fileID,'#Stacking option\n');
fprintf(fileID, sprintf('Stacking_Frequency     =%8.4f\n',Stacking_Frequency));
fprintf(fileID, sprintf('Source_average_num_per_unithour=%d\n',Source_average_num_per_unithour));
fprintf(fileID, sprintf('C3_truncate_alpha =%8.4f\n',C3_truncate_alpha));
fprintf(fileID, sprintf('C3_max_time_lag =%8.4f\n',C3_max_time_lag));
fprintf(fileID,'#Output option\n');
fprintf(fileID, sprintf('IsSaveRawTrace         =%s\n',IsSaveRawTrace));
fprintf(fileID, sprintf('IsSaveCC               =%s\n',IsSaveCC));
fclose(fileID);

