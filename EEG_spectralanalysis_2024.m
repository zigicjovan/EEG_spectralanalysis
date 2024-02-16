%{
2024/02/11, Jovan Zigic, Dominican University New York

This script generates spectral (frequency) analysis results for an EEG
experiment.

Required input parameters: 
Fs := EEG sampling rate (Hz) a.k.a. number of samples per second
numberofchannels := number of EEG channels
subjects := number of subjects in the experiment
numberofdatasets := number of datasets being analyzed
dataset* := the file name of each data set for each subject in the experiment
instructions := chosen samples of raw data to use in analysis

Work in progress or existing issues (ranked by priority):
1. Mean output bar chart must be adapted to only include cleaned data
2. Figure legends should be adapted for selected channels
3. Can bad data be salvaged?
4. Loop and output should be adapted for various number of data sets
5. What is the theoretical criteria for usable data (check spectrum decay)?
%}

tic % start timer for entire script (runtime displayed at the end of the script)

clear;
close all;

% Experiment parameters
Fs = 300; % sampling rate (Hz) a.k.a. number of samples per second
numberofchannels_raw = 7; % number of EEG channels
subjects = 20; % number of subjects in experiment

% file names to analyze
numberofdatasets = 2;
dataset1 = 'pre-test.csv';
dataset2 = 'post-test.csv'; 

% cleaning instructions
cleaning_instructions  = readmatrix(['EEG_cleandata_2024' ...
    '.xlsx']);
cleaning_instructions = cleaning_instructions( : , 2:end); % remove subject number from matrix

% frequency band ranges and names
bands = [ 0.5,4 ; 4.5,8 ; 8.5,13 ; 13.5,30 ; 30.5,100 ]; 
numberofbands = size(bands,1);
bandnames = categorical([ "Delta", "Theta", "Alpha", "Beta", "Gamma" ]);

% output file with frequency band distribution for each subject 
% (row := subject ; column := {dataset1(d,t,a,b,g)}, {dataset2(d,t,a,b,g)}, ...)
normalized_percent_bandpower = zeros(subjects,numberofbands*numberofdatasets);

for sub = 1:subjects

    %%% Read voltage-time data and convert to frequency data %%%
    
    % Prepare data sets for Pre and Post
    file1 = [ pwd '\Participant_Data_2024\subject ' num2str(sub) ' - ' dataset1 '']; % file path
    file2 = [ pwd '\Participant_Data_2024\subject ' num2str(sub) ' - ' dataset2 ''];
    datatable1  = readtable(file1,'VariableNamingRule','preserve'); % read csv file
    datatable2  = readtable(file2,'VariableNamingRule','preserve'); 
    t1_raw = table2array(datatable1(:,1)); % convert to array format
    t2_raw = table2array(datatable2(:,1));
    voltagedata1_raw = table2array(datatable1(:,2:(numberofchannels_raw + 1))); % isolate channels from data file
    voltagedata2_raw = table2array(datatable2(:,2:(numberofchannels_raw + 1))); 

    %{
    % Prepare data sets for Tx1 and Tx2
    fileA = [ pwd '\Participant Data\Subject ' num2str(sub) '\' datasetA '']; 
    fileB = [ pwd '\Participant Data\Subject ' num2str(sub) '\' datasetB ''];
    datatableA  = readtable(fileA,'VariableNamingRule','preserve'); 
    datatableB  = readtable(fileB,'VariableNamingRule','preserve'); 
    tA_raw = table2array(datatableA(:,1)); 
    tB_raw = table2array(datatableB(:,1));
    voltagedataA_raw = table2array(datatableA(:,2:(numberofchannels_raw + 1))); 
    voltagedataB_raw = table2array(datatableB(:,2:(numberofchannels_raw + 1))); 
    %}

    %%% clean data %%%
    try % choose which segments of the data to use
        [ t1 , voltagedata1 ] = clean_data( t1_raw , voltagedata1_raw , ...
            cleaning_instructions( sub, ((numberofchannels_raw+2)*0 + 1):(numberofchannels_raw+2)*1) );
        numberofchannels_1 = size( voltagedata1 , 2 );
        voltagedata1(1,1); % will catch error if empty => use entire data set when no cleaning is specified
    catch % use raw data without cleaning
        t1 = t1_raw;
        voltagedata1 = voltagedata1_raw;
        numberofchannels_1 = size( voltagedata1 , 2 );
    end

    try
        [ t2 , voltagedata2 ] = clean_data( t2_raw , voltagedata2_raw , ...
            cleaning_instructions( sub, ((numberofchannels_raw+2)*1 + 1):end ) );
        numberofchannels_2 = size( voltagedata2 , 2 );
        voltagedata2(1,1); 
    catch
        t2 = t2_raw;
        voltagedata2 = voltagedata2_raw;
        numberofchannels_2 = size( voltagedata2 , 2 );
    end

    %{
    try
        [ tA , voltagedataA ] = clean_data( tA_raw , voltagedataA_raw , ...
            cleaning_instructions( sub, ((numberofchannels_raw+2)*1 + 1):(numberofchannels_raw+2)*2) );
        numberofchannels_A = size( voltagedataA , 2 );
        voltagedataA(1,1); 
    catch
        tA = tA_raw;
        voltagedataA = voltagedataA_raw;
        numberofchannels_A = size( voltagedataA , 2 );
    end
    
    try
        [ tB , voltagedataB ] = clean_data( tB_raw , voltagedataB_raw , ...
            cleaning_instructions( sub, ((numberofchannels_raw+2)*2 + 1):(numberofchannels_raw+2)*3) );
        numberofchannels_B = size( voltagedataB , 2 );
        voltagedataB(1,1); 
    catch
        tB = tB_raw;
        voltagedataB = voltagedataB_raw;
        numberofchannels_B = size( voltagedataB , 2 );
    end
    %}

    %
    % inspect data visually
    figure(1)
    for j = 1:numberofchannels_1
        plot(t1,voltagedata1(:,j))
        hold on
    end
    xlabel('Time (seconds)');
    ylabel('Voltage (microvolts)');
    title(['EEG Readings Pre Test for Subject ' num2str(sub) ]);
    legend({'LE','F4','C4','P4','P3','C3','F3'},'Location','southoutside', 'NumColumns', 7)
    hold off

    figure(2)
    for j = 1:numberofchannels_2
        plot(t2,voltagedata2(:,j))
        hold on
    end
    xlabel('Time (seconds)');
    ylabel('Voltage (microvolts)');
    title(['EEG Readings Post Test for Subject ' num2str(sub) ]);
    legend({'LE','F4','C4','P4','P3','C3','F3'},'Location','southoutside', 'NumColumns', 7)
    hold off

    %{
    figure(3)
    for j = 1:numberofchannels_A
        plot(tA,voltagedataA(:,j))
        hold on
    end
    xlabel('Time (seconds)');
    ylabel('Voltage (microvolts)');
    title(['EEG Readings during PT Treatment 1 for Subject ' num2str(sub) ]);
    legend({'LE','F4','C4','P4','P3','C3','F3'},'Location','southoutside', 'NumColumns', 7)
    hold off

    figure(4)
    for j = 1:numberofchannels_B
        plot(tB,voltagedataB(:,j))
        hold on
    end
    xlabel('Time (seconds)');
    ylabel('Voltage (microvolts)');
    title(['EEG Readings during PT Treatment 2 for Subject ' num2str(sub) ]);
    legend({'LE','F4','C4','P4','P3','C3','F3'},'Location','southoutside', 'NumColumns', 7)
    hold off
    %}

    % transform data from time (voltage) domain to frequency (fourier) domain
    freqdata1 = real(fft(voltagedata1));
    freqdata2 = real(fft(voltagedata2));
    %freqdataA = real(fft(voltagedataA));
    %freqdataB = real(fft(voltagedataB));
    
    %%% Perform spectral analysis of signal %%%
    
    % initialize variables for frequency bands
    percent_bandpower1 = zeros(1,numberofbands);
    pBand1 = zeros(numberofbands,numberofchannels_1);
    pTot1 = zeros(numberofbands,numberofchannels_1);
    percent_bandpower2 = zeros(1,numberofbands);
    pBand2 = zeros(numberofbands,numberofchannels_2);
    pTot2 = zeros(numberofbands,numberofchannels_2);
    %{
    percent_bandpowerA = zeros(1,numberofbands);
    pBandA = zeros(numberofbands,numberofchannels_A);
    pTotA = zeros(numberofbands,numberofchannels_A);
    percent_bandpowerB = zeros(1,numberofbands);
    pBandB = zeros(numberofbands,numberofchannels_B);
    pTotB = zeros(numberofbands,numberofchannels_B);
    %}  
    for i = 1:numberofbands % perform for each frequency band
        x1 = freqdata1;
        [Pxx1,F1] = pwelch(x1,[],[],[],Fs); % Welch's power spectral density estimate (8 segments with 50% overlap)
        %[Pxx1,F1] = periodogram(x1,rectwin(length(x1)),length(x1),Fs); % power spectral density estimate
        pBand1(i,:) = bandpower(Pxx1,F1,bands(i,:),'psd'); % power of signal in frequency band
        pTot1(i,:) = bandpower(Pxx1,F1,'psd'); % power of signal
        percent_bandpower1(1,i) = 100*(pBand1(i,:)/pTot1(i,:)); % proportional power of frequency band

        x2 = freqdata2;
        [Pxx2,F2] = pwelch(x2,[],[],[],Fs); 
        %[Pxx2,F2] = periodogram(x2,rectwin(length(x2)),length(x2),Fs); 
        pBand2(i,:) = bandpower(Pxx2,F2,bands(i,:),'psd'); 
        pTot2(i,:) = bandpower(Pxx2,F2,'psd'); 
        percent_bandpower2(1,i) = 100*(pBand2(i,:)/pTot2(i,:)); 
        %{
        xA = freqdataA;
        [PxxA,FA] = pwelch(xA,[],[],[],Fs); 
        %[PxxA,FA] = periodogram(xA,rectwin(length(xA)),length(xA),Fs); 
        pBandA(i,:) = bandpower(PxxA,FA,bands(i,:),'psd');
        pTotA(i,:) = bandpower(PxxA,FA,'psd'); 
        percent_bandpowerA(1,i) = 100*(pBandA(i,:)/pTotA(i,:)); 

        xB = freqdataB;
        [PxxB,FB] = pwelch(xB,[],[],[],Fs);
        %[PxxB,FB] = periodogram(xB,rectwin(length(xB)),length(xB),Fs); 
        pBandB(i,:) = bandpower(PxxB,FB,bands(i,:),'psd'); 
        pTotB(i,:) = bandpower(PxxB,FB,'psd'); 
        percent_bandpowerB(1,i) = 100*(pBandB(i,:)/pTotB(i,:)); 
        %}
    end
    
    % normalized proportion of frequency band
    normalized_percent_bandpower( sub, (numberofbands*0)+1 : (numberofbands*1) ) = percent_bandpower1/sum(percent_bandpower1);
    npb1 = normalized_percent_bandpower( sub, (numberofbands*0)+1 : (numberofbands*1) );
    normalized_percent_bandpower( sub, (numberofbands*1)+1 : (numberofbands*2) ) = percent_bandpower2/sum(percent_bandpower2);
    npb2 = normalized_percent_bandpower( sub, (numberofbands*1)+1 : (numberofbands*2) );
    %{
    normalized_percent_bandpower( sub, (numberofbands*2)+1 : (numberofbands*3) ) = percent_bandpowerA/sum(percent_bandpowerA);
    npb3 = normalized_percent_bandpower( sub, (numberofbands*2)+1 : (numberofbands*3) );
    normalized_percent_bandpower( sub, (numberofbands*3)+1 : (numberofbands*4) ) = percent_bandpowerB/sum(percent_bandpowerB);
    npb4 = normalized_percent_bandpower( sub, (numberofbands*3)+1 : (numberofbands*4) );
    %}
    subject_bar_chart_values = [ npb1 ; npb2]; 

    % visualize distribution of signal band power a.k.a. power spectral density
    %
    figure(5)
    bar(bandnames,subject_bar_chart_values)
    xlabel('EEG Frequency Bands');
    ylabel('Power Spectral Density');
    title(['EEG Signal Relative Power Spectral Density for Subject ' num2str(sub) ]);
    legend({'Pre Test','Post Test'},'Location','southoutside', 'NumColumns', 7)
    %}

    %%% save files to current directory %%%
    %
    filename1 = [ 'EEGreadings_pre_subject' num2str(sub) ];
    saveas(figure(1),filename1,'png');
    filename2 = [ 'EEGreadings_post_subject' num2str(sub) ];
    saveas(figure(2),filename2,'png');
    %{
    filename3 = [ 'EEGreadings_tx1_subject' num2str(sub) ];
    saveas(figure(3),filename3,'png');
    filename4 = [ 'EEGreadings_tx2_subject' num2str(sub) ];
    saveas(figure(4),filename4,'png');
    %}
    filename5 = [ 'relativePSD_subject' num2str(sub) '' ];
    saveas(figure(5),filename5,'png');
    %}
end

%%% Post-processing of all experimental data %%%

% Visualize average distribution of power spectral density across subjects
mean_npb = zeros( 1, numberofbands*numberofdatasets );
for j = 1:numberofbands*numberofdatasets
    mean_npb(1,j) = mean(normalized_percent_bandpower( : , j ));    
end
mean_bar_chart_values = [ mean_npb(1,(numberofbands*0)+1 : (numberofbands*1)) ; 
                        mean_npb(1,(numberofbands*1)+1 : (numberofbands*2))  ]; 
figure(6)
bar(bandnames,mean_bar_chart_values)
xlabel('EEG Frequency Bands');
ylabel('Power Spectral Density');
title(['Average Relative PSD for ' num2str(subjects) ' subjects' ]);
legend({'Pre Test','Post Test'},'Location','southoutside', 'NumColumns', 7)
filename6 = 'mean_relativePSD_allsubjects' ;
saveas(figure(6),filename6,'png');

% Save data to CSV file
ColumnNames = ["Subject","PreDelta","PreTheta",'PreAlpha','PreBeta','PreGamma',...
    'PostDelta','PostTheta','PostAlpha','PostBeta','PostGamma'];
SubjectNumber = zeros(subjects,1);
for j = 1:subjects
    SubjectNumber(j,1) = j;
end
subjectanddata = [ SubjectNumber , normalized_percent_bandpower];
writematrix(ColumnNames,'relative_PSD_allsubjects.csv');
writematrix(subjectanddata,'relative_PSD_allsubjects.csv','WriteMode','append');

% stop and display timer for entire script
runtime = toc; 
disp(['This took ' num2str(runtime) ' seconds to compute!'])

function [ t , data ] = clean_data( t_raw , data_raw , instructions)

    % this function cleans the data from artifacts or uses only desired segments of the raw data set
    
    % 1. time segment
    ta = instructions(1,1);
    tb = instructions(1,2);
    t_segment = [ ta , tb ];

    t_start = max(find( t_raw <= t_segment(1) ));
    if isempty(t_start)
        t_start = 1;
    end
    t_end = max(find( t_raw <= t_segment(2) ));
    t = t_raw( t_start:t_end , 1);
    data = data_raw( t_start:t_end , : );

    % 2. remove channels (if necessary) 
    channel_segment = instructions(1,3:end);
    channel_segment = rmmissing(channel_segment);
    data = data( : , channel_segment );

end