%%% Routine for Electromiographic Data and Postural Automatism.%%%
%%% By André Ivaniski Mello & Lucas Liz Alves %%%
%%% andreivaniskimello@gmail.com, lucaslizalves@hotmail.com %%%


%%%INPUT: emg .txt file (in format of two columns. First column for time.
%%%Second column for EMG activity. One muscle channel below the other. Channel is equal to EMG muscle, because ONE muscle is collected in ONE EMG Channel%%%

%%%% OUTPUT: BASELINE,APA and CPA EMG activity of Postural Muscles. The Output is an .csv (Excel) file. %%%%

%ATTENTION!! Before run the routine, some variables need to be adjusted:

% Subject name
% Condition name
% File directory with the Raw EMG data

% Butterworth filters parameters (cut off frequencies, order) (default is
% PassBand for Agonist Muscle of 40 to 450 Hz 2º order; a Low Pass Band of
% 40 Hz, 2º order for Postural Muscles)

% Sample frequency (fs) (default is 2000 Hz)
% RMS window lenght desired (default is 100 ms)
% Threshold for t0 determination from Agonist RMS activity (default is
% Mean of Basal Activity*2)

% Columns Numbers of Postural Muscles (default is from 5 to 7)

% Baseline, APA and CPA period of time in relation to t0 (Default is
% Baseline: -500 to -400 ms; APA: -100 to 0 ms; CPA: 0 to +100 ms)

% After running the Routine for the first time, visually check if the t0
% values extracted by threshold method are good. If not, if they do not
% correspond well to the Agonist RMS curve, manually fill the "Manual t0
% Matrix" in the t0 SECTION.

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%INFORMATION FOR DATA EXPORTATION SECTION

%Insert Output File Path here
File_Name =['Camila_ABDUNP_TESTE']; 
Output_File_Path = ['C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto APA\Data Out'];   

%%
%%%% RAW DATA FILES DIRECTORIES SECTION %%%%

%%% Load EMG data of TWO EMGs files (each one with 08 channels) and Separate each Muscle Channel

%Insert Raw Data Files directories here:

file_1 = 'C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto APA\Dados APA EMG Solo\TXT\Camila\CAABDUNP1.PRN';    %EMG 01 with the Agonist in First Channel

file_2 = 'C:\Users\andre\Documents\Mestrado\Projetos Pesquisa\Projeto APA\Dados APA EMG Solo\TXT\Camila\CAABDUNP2.PRN';    %EMG 02

%% LOAD EMG FILE 01 SECTION

fid = fopen(file_1,'r');
data_full_1 = textscan(fid, '%s%s', 'Delimiter',' ', 'headerLines', 13);    %Read all data as string
fclose(fid);

fid = fopen(file_1,'r');
data_time_1 = textscan(fid, '%f%f', 'Delimiter',' ', 'headerLines', 13);    %Read only the numerical values (n x 2) of the First Channel
fclose(fid);

time_1 = data_full_1{1};            %Returns Time vector (n x 1) of all channels. One channel after the other in one column.
emg_1 = data_full_1{2};             %Returns EMG vector (n x 1) of all channels. One channel after the other in one column.

Time_Vector_1 = data_time_1{1};     %Return Time vector of only the first channel.
size_1 = size(Time_Vector_1);       %Return the size (rows) of Time vector of the first channel. This value it will be used the recort the big time and emg vector per channels.
size_1 = size_1(1);


%NAME OF THE CHANNELS SECTION

fid = fopen(file_1,'r');
data_full_1 = textscan(fid, '%s%s', 'Delimiter',' ', 'headerLines', 9); %read all data as string
fclose(fid);
name_1 = data_full_1{1};        %Return the first column of the full file (with headers inclusive)


Channel_1_Name = name_1(1);   %Return the Name of the Channel 1. This is out of loop, because it do not have to be multiply be the size of Time vector.


%This is a loop to extract the name of each EMG Channel. Channels_Name is
%the matrix with the name of each Channel

for j = 1:(((length(emg_1))/size_1))-1
    Channels_Name_1(j) = name_1((size_1*j)+(5+4*(j-1)));                %Return the Name of the Channel 2. This summation corrects for the offset between Time Vector Size and Full File with Headers.
    Number_of_Channels_1 = j+1;               %Return the Number of EMG Channels just to check
end

Channels_Name_1 = [Channel_1_Name Channels_Name_1];

Channel_1_Name = Channels_Name_1(1);          %% ATENTION: ADJUST THIS PART MANUALLY FOR THE NUMBERS OF CHANNELS. Made a Loop for this function is not reccomended (https://www.mathworks.com/matlabcentral/answers/304528-tutorial-why-variables-should-not-be-named-dynamically-eval#answer_236124)
Channel_2_Name = Channels_Name_1(2); 
Channel_3_Name = Channels_Name_1(3);
Channel_4_Name = Channels_Name_1(4);
Channel_5_Name = Channels_Name_1(5);
Channel_6_Name = Channels_Name_1(6);
Channel_7_Name = Channels_Name_1(7);
Channel_8_Name = Channels_Name_1(8);


%SEPARATE EACH MUSCLE SECTION
%This part is spliting the big row vector array in individual channels.

Channel_1_Time = time_1(1:size_1,1);
Channel_1_Emg = emg_1(1:size_1,1);

%for i = 1:(Number_of_Channels-1)
%    Channels_EMG(i) = emg(((i*size)+(5+4*(i-1)):((i+1)*size)));
%end

Channel_2_Emg = emg_1(((1*size_1)+5):((2*size_1)+(4*1)));               %size+5: because each channel has a Headerlines of 4 rows. Size = Channel_n last line. Size+5 = Channel_n+1 first line.
Channel_3_Emg = emg_1(((2*size_1)+(5+(4*1))):((3*size_1)+(4*2)));
Channel_4_Emg = emg_1(((3*size_1)+(5+(4*2))):((4*size_1)+(4*3)));
Channel_5_Emg = emg_1(((4*size_1)+(5+(4*3))):((5*size_1)+(4*4)));
Channel_6_Emg = emg_1(((5*size_1)+(5+(4*4))):((6*size_1)+(4*5)));
Channel_7_Emg = emg_1(((6*size_1)+(5+(4*5))):((7*size_1)+(4*6)));
Channel_8_Emg = emg_1(((7*size_1)+(5+(4*6))):((8*size_1)+(4*7)));


%This part is converting the string array type into numeric array type.

Channel_1_Emg = str2double(Channel_1_Emg);
Channel_2_Emg = str2double(Channel_2_Emg);
Channel_3_Emg = str2double(Channel_3_Emg);
Channel_4_Emg = str2double(Channel_4_Emg);
Channel_5_Emg = str2double(Channel_5_Emg);
Channel_6_Emg = str2double(Channel_6_Emg);
Channel_7_Emg = str2double(Channel_7_Emg);
Channel_8_Emg = str2double(Channel_8_Emg);

All_Channels_Emg_1 = [Channel_1_Emg Channel_2_Emg Channel_3_Emg Channel_4_Emg Channel_5_Emg Channel_6_Emg Channel_7_Emg Channel_8_Emg];

%% LOAD EMG FILE 02 SECTION

fid = fopen(file_2,'r');
data_full_2 = textscan(fid, '%s%s', 'Delimiter',' ', 'headerLines', 13); %read all data as string
fclose(fid);

fid = fopen(file_2,'r');
data_time_2 = textscan(fid, '%f%f', 'Delimiter',' ', 'headerLines', 13);
fclose(fid);

time_2 = data_full_2{1}; %Returns Time vector (n x 1) of all channels. One channel after the other in one column.
emg_2 = data_full_2{2};  %Returns EMG vector (n x 1) of all channels. One channel after the other in one column.

Time_Vector_2 = data_time_2{1}; %Return Time vector of only the first channel.
size_2 = size(Time_Vector_2); %Return the size (rows) of Time vector of the first channel. This value it will be used the recort the big time and emg vector per channels.
size_2 = size_2(1);


%NAME OF THE CHANNELS SECTION

fid = fopen(file_2,'r');
data_full_2 = textscan(fid, '%s%s', 'Delimiter',' ', 'headerLines', 9); %read all data as string
fclose(fid);
name_2 = data_full_2{1};        %Return the first column of the full file (with headers inclusive)


Channel_9_Name = name_2(1);   %Return the Name of the Channel 1. This is out of loop, because it do not have to be multiply be the size of Time vector.


%This is a loop to extract the name of each EMG Channel. Channels_Name is
%the matrix with the name of each Channel

for j = 1:(((length(emg_2))/size_2))-1
    Channels_Name_2(j) = name_2((size_2*j)+(5+4*(j-1)));        %Return the Name of the Channel 2. This summation corrects for the offset between Time Vector Size and Full File with Headers.
    Number_of_Channels_2 = j+1;                                 %Return the Number of EMG Channels just to check
end

Channels_Name_2 = [Channel_9_Name Channels_Name_2];

Channel_9_Name = Channels_Name_2(1);          %% ATENTION: ADJUST THIS PART MANUALLY FOR THE NUMBERS OF CHANNELS. Made a Loop for this function is not reccomended (https://www.mathworks.com/matlabcentral/answers/304528-tutorial-why-variables-should-not-be-named-dynamically-eval#answer_236124)
Channel_10_Name = Channels_Name_2(2); 
Channel_11_Name = Channels_Name_2(3);
Channel_12_Name = Channels_Name_2(4);
Channel_13_Name = Channels_Name_2(5);
Channel_14_Name = Channels_Name_2(6);
Channel_15_Name = Channels_Name_2(7);
Channel_16_Name = Channels_Name_2(8);


%SEPARATE EACH MUSCLE SECTION
%This part is spliting the big row vector array in individual channels.

Channel_9_Time = time_2(1:size_2,1);
Channel_9_Emg = emg_2(1:size_2,1);

%for i = 1:(Number_of_Channels-1)
%    Channels_EMG(i) = emg(((i*size)+(5+4*(i-1)):((i+1)*size)));
%end

Channel_10_Emg = emg_2(((1*size_2)+5):((2*size_2)+(4*1)));               %size+5: because each channel has a Headerlines of 4 rows. Size = Channel_n last line. Size+5 = Channel_n+1 first line.
Channel_11_Emg = emg_2(((2*size_2)+(5+(4*1))):((3*size_2)+(4*2)));
Channel_12_Emg = emg_2(((3*size_2)+(5+(4*2))):((4*size_2)+(4*3)));
Channel_13_Emg = emg_2(((4*size_2)+(5+(4*3))):((5*size_2)+(4*4)));
Channel_14_Emg = emg_2(((5*size_2)+(5+(4*4))):((6*size_2)+(4*5)));
Channel_15_Emg = emg_2(((6*size_2)+(5+(4*5))):((7*size_2)+(4*6)));
Channel_16_Emg = emg_2(((7*size_2)+(5+(4*6))):((8*size_2)+(4*7)));


%This part is converting the string array type into numeric array type.

Channel_9_Emg = str2double(Channel_9_Emg);
Channel_10_Emg = str2double(Channel_10_Emg);
Channel_11_Emg = str2double(Channel_11_Emg);
Channel_12_Emg = str2double(Channel_12_Emg);
Channel_13_Emg = str2double(Channel_13_Emg);
Channel_14_Emg = str2double(Channel_14_Emg);
Channel_15_Emg = str2double(Channel_15_Emg);
Channel_16_Emg = str2double(Channel_16_Emg);

All_Channels_Emg_2 = [Channel_9_Emg Channel_10_Emg Channel_11_Emg Channel_12_Emg Channel_13_Emg Channel_14_Emg Channel_15_Emg Channel_16_Emg];


figure(1)
plot(All_Channels_Emg_1(:,3),'m')
hold on
plot(All_Channels_Emg_1(:,8),'g')
hold on
plot(All_Channels_Emg_2(:,8),'r')
title ('EMG');
xlabel('time')
ylabel('Emg')
legend('DMD', 'Sync 01', 'Sync 02')

%%
%%%% SYNC ADJUSTMENT SECTION %%%%

Sync_Basal_1 = abs(mean(Channel_8_Emg(1:500)));
    
    if Sync_Basal_1 < 1
    Sync_1 = find(Channel_8_Emg > 1000*Sync_Basal_1);
    elseif Sync_Basal_1 >= 1
    Sync_1 = find(Channel_8_Emg > 100*Sync_Basal_1);

    end
Sync_1 = Sync_1(1);

All_Channels_Emg_1 = All_Channels_Emg_1(Sync_1:end,:);

Sync_Basal_2 = abs(mean(Channel_16_Emg(1:500)));
 
if Sync_Basal_2 < 1
    Sync_2 = find(Channel_16_Emg > 1000*Sync_Basal_2);
 elseif Sync_Basal_2 >= 1
    Sync_2 = find(Channel_16_Emg > 100*Sync_Basal_2);
 end
 
Sync_2 = find(Channel_16_Emg > 100*Sync_Basal_2);
Sync_2 = Sync_2(1);

All_Channels_Emg_2 = All_Channels_Emg_2(Sync_2:end,:);

Length_Emg = min(length(All_Channels_Emg_1), length(All_Channels_Emg_2));

All_Channels_Emg_1 = All_Channels_Emg_1(1:Length_Emg,:);
All_Channels_Emg_2 = All_Channels_Emg_2(1:Length_Emg,:);

fs = 2000;
Time_After_Sync = size(All_Channels_Emg_1);
Time_After_Sync = Time_After_Sync(1);
Time_After_Sync =  0:(1/fs):((Time_After_Sync-1)/fs);
Time_After_Sync = Time_After_Sync.';

%% Create ONE Matriz for all 16 Muscles of the 02 EMGs files


All_Channels_Emg = [All_Channels_Emg_1 All_Channels_Emg_2];
Number_of_Channels = Number_of_Channels_1 + Number_of_Channels_2;

Channels_Name = [Channels_Name_1 Channels_Name_2];
%%

%%%% DATA FILTERING SECTION %%%%

%%%% Remove Offset of the EMG signal %%%%

%The loop is to Remove Offset of all
%muscles (One muscle per Column)
for j=1:16
    All_Channels_Emg_Offset(:,j) = detrend(All_Channels_Emg(:,j),0);
end


%%%% Filter Butterworth Passband of the EMG signal %%%%

%The loop is to Passband Filter all muscles (One muscle per Column)

%ATENTION: Adjust Filter Order (order), Low (fcutlow) and High (fcuthigh) Cutoff Frequencies, and Sample Frequency
%(fs)

fs = 2000;
order = 2 ;
fcutlow  = 40;
fcuthigh = 400;
[b,a] = butter(order,[fcutlow,fcuthigh]/(fs/2), 'bandpass');
  
for j=1:Number_of_Channels
     All_Channels_Emg_Butter(:,j) = filtfilt(b,a,All_Channels_Emg_Offset(:,j));
end

% Here are TWO STOP BAND BUTTERWORTH filters to remove 60 and 120 Hz
% frequencies from all the EMG signals

fcutlow_stopband  = 59;
fcuthigh_stopband = 61;
[b_stop,a_stop] = butter(order,[fcutlow_stopband,fcuthigh_stopband]/(fs/2), 'stop');

for j=1:Number_of_Channels
     All_Channels_Emg_Butter(:,j) = filtfilt(b_stop,a_stop,All_Channels_Emg_Butter(:,j));
end

fcutlow_stopband  = 119;
fcuthigh_stopband = 121;
[b_stop,a_stop] = butter(order,[fcutlow_stopband,fcuthigh_stopband]/(fs/2), 'stop');

for j=1:Number_of_Channels
     All_Channels_Emg_Butter(:,j) = filtfilt(b_stop,a_stop,All_Channels_Emg_Butter(:,j));
end


%%
%%%%% RMS of AGONIST SECTION %%%%%

Agonist =  All_Channels_Emg_Butter(:,3);        % Line to Label the agonist Muscle of the movement

%Insert the RMS window length desired in ms below
windowlength_in_ms = 100;                                       %in ms
windowlength_in_lines = (windowlength_in_ms/1000)*fs;           %convert windowlength into number of lines accordingly to the sample frequency

Agonist = Agonist.^2;               %Square of agonist to calculate RMS inside the loop below

for j = 1:(length(Agonist)/windowlength_in_lines)
	Agonist_RMS(j) = sqrt(mean(Agonist(1+windowlength_in_lines*(j-1):windowlength_in_lines*(j)-1)));
    
end

Agonist_RMS = Agonist_RMS.';        %Transpose matrix of One Line x n columns into a matrix of n Lines x One Column.


figure (2)
title ('EMG Activity');
xlabel('Time (s)');
ylabel ('Voltage (mV)');
hold on
plot (Agonist_RMS,'g');
hold on
legend('Butter')
%%
%%%% t0 SECTION %%%%
    
    %%% Determination (repetitions start) %%%%
  
    %OBSERVATION ABOUT t0 DETERMINATION: Run the routine for the first
    %time. Check visually if the t0 matrix extracted by threshold method
    %corresponds to the real t0. If not, manually set the t0 values in the
    %matrix below.

%Set the Basal Activity period here!!

Basal_Start = 2;
Basal_End = 40;

Agonist_Basal_Mean = mean(Agonist_RMS(Basal_Start:Basal_End));
Agonist_Basal_StdDev = std(Agonist_RMS(Basal_Start:Basal_End));

% Threshold detector of Agonist RMS activity. 
% The literature says to use a threshold of Mean + 2 SD; but this threshold
% was too low, not sensitive. So I set a threshold of 2*Mean. Mean and SD
% of the Basal Activity calculated above.

t_threshold = find(Agonist_RMS > (Agonist_Basal_Mean*4));       


for j = 1:(length(t_threshold)-1)
    if t_threshold(j+1) - t_threshold(j) > 1
        tend(j) = t_threshold(j);
        t0(j) = t_threshold(j+1);
   end 
end

tend = tend.';      % End of Repetitions
tend(tend==0)=[];   % Delete zero elements from matrix

t0 = t0.';          
t0(t0==0)=[];
t0 = [t_threshold(1); t0];  % Start (t0) of Repetitions
t0(t0<Basal_End)=[];

%"Manual t0 Matrix". UNBLOCK THIS SECTION IF NECESSARY!!

%Set here the t0 values manually after visual
%inspection of the Agonist RMS curve, if the threshold t0 is not
%corresponding well to the real values. 

t0_slow = [54; 83; 105];
t0_fast = [173; 185; 196];
t0_load = [288; 313; 335];
t0 = [t0_slow; t0_fast; t0_load];
% End of "Manual t0 Matrix"

t0_in_s = t0 * windowlength_in_lines * (1/fs);
t0_in_lines = t0_in_s * fs;
[nlt0 nct0] = size(t0_in_lines);                        %nlt0 is equal to the number of repetitions.

%%

%%%%   POSTURAL MUSCLES FILTER SECTION %%%%

%%%Adjust the Columns of Postural Muscles

First_Postural_Muscle_Column_Emg_1 = 5;
Last_Postural_Muscle_Column_Emg_1 = 7;

First_Postural_Muscle_Column_Emg_2 = 9;
Last_Postural_Muscle_Column_Emg_2 = 15;

Postural_Muscles_Name_1 = Channels_Name(:,First_Postural_Muscle_Column_Emg_1:Last_Postural_Muscle_Column_Emg_1); 
Postural_Muscles_Name_2 = Channels_Name(:,First_Postural_Muscle_Column_Emg_2:Last_Postural_Muscle_Column_Emg_2);

Postural_Muscles_Name = [Postural_Muscles_Name_1 Postural_Muscles_Name_2];

Postural_Muscles_Emg_1 = All_Channels_Emg_Offset(:,First_Postural_Muscle_Column_Emg_1:Last_Postural_Muscle_Column_Emg_1);
Postural_Muscles_Emg_2 = All_Channels_Emg_Offset(:,First_Postural_Muscle_Column_Emg_2:Last_Postural_Muscle_Column_Emg_2);

Postural_Muscles_Emg_Integral = [Postural_Muscles_Emg_1 Postural_Muscles_Emg_2];

[nl_Postural nc_Postural] = size(Postural_Muscles_Emg_Integral);

%EMG of Postural Muscles rectification (absolute values)
Postural_Muscles_Emg_Rectified = abs(Postural_Muscles_Emg_Integral);

%EMG of Postural Muscles Low Pass Butterworth Filter
fs = 2000;
order = 2 ;
fcuthigh = 40;
[b,a] = butter(order,fcuthigh/(fs/2), 'low');
  
for j=1:nc_Postural 
     Postural_Muscles_Emg_Butter(:,j) = filtfilt(b,a,Postural_Muscles_Emg_Rectified(:,j));
end


fcutlow_stopband  = 59;
fcuthigh_stopband = 61;
[b_stop,a_stop] = butter(order,[fcutlow_stopband,fcuthigh_stopband]/(fs/2), 'stop');

for j=1:nc_Postural
     Postural_Muscles_Emg_Butter(:,j) = filtfilt(b_stop,a_stop,Postural_Muscles_Emg_Butter(:,j));
end


fcutlow_stopband  = 119;
fcuthigh_stopband = 121;
[b_stop,a_stop] = butter(order,[fcutlow_stopband,fcuthigh_stopband]/(fs/2), 'stop');

for j=1:nc_Postural 
     Postural_Muscles_Emg_Butter(:,j) = filtfilt(b_stop,a_stop,Postural_Muscles_Emg_Butter(:,j));
end

%%

%%%% APA and CPA ANALYSIS SECTION %%%%

% Postural Automatisms Windows from: Kanekar & Aruin, 2015

%%% APA = -100 ms to +50
%%% CPA = +50 to +200 ms
%%% Baseline Activity = -500 ms to -400ms

APA_in_ms_Start = 100;
APA_in_ms_End = 50;

CPA_in_ms_Start = 50;
CPA_in_ms_End = 200;

Baseline_in_ms_Start = 500;
Baseline_in_ms_End = 400;


APA_in_lines_Start = (APA_in_ms_Start/1000) * fs;
APA_in_lines_End = (APA_in_ms_End/1000) * fs;

CPA_in_lines_Start = (CPA_in_ms_Start/1000) *  fs;
CPA_in_lines_End = (CPA_in_ms_End/1000) *  fs;

Baseline_in_lines_Start = (Baseline_in_ms_Start/1000) * fs;
Baseline_in_lines_End = (Baseline_in_ms_End/1000) * fs;


% INTEGRAL CALCULATION OF EACH POSTURAL AUTOMATISM TIME WINDOW FOR EACH
% POSTURAL MUSCLE

windowlength_postural_in_ms = 10;                                       %in ms
windowlength_postural_in_lines = (windowlength_postural_in_ms/1000)*fs;


for i = 1:nlt0
Postural_Muscles_Emg_Baseline_Window = Postural_Muscles_Emg_Butter(t0_in_lines(i) - Baseline_in_lines_Start : t0_in_lines(i) - Baseline_in_lines_End, :);
Postural_Muscles_Emg_APA_Window = Postural_Muscles_Emg_Butter(t0_in_lines(i) - APA_in_lines_Start : t0_in_lines(i) + APA_in_lines_End, :);
Postural_Muscles_Emg_CPA_Window = Postural_Muscles_Emg_Butter(t0_in_lines(i) + CPA_in_lines_Start : t0_in_lines(i) + CPA_in_lines_End, :);

	% Loop to calculate INTEGRAL of Baseline, APA, CPA 
    for j=1:nc_Postural 
    Postural_Muscles_Emg_Baseline_Integral_Vector(:,j) = cumtrapz(Postural_Muscles_Emg_Baseline_Window(:,j));
    Postural_Muscles_Emg_APA_Integral_Vector(:,j) = cumtrapz(Postural_Muscles_Emg_APA_Window(:,j));
    Postural_Muscles_Emg_CPA_Integral_Vector(:,j) = cumtrapz(Postural_Muscles_Emg_CPA_Window(:,j));
    end
    
    %Here the last value of the INTEGRAL vector is being extracted for each postural muscle
    Postural_Muscles_Emg_Baseline_Integral(i,:) = Postural_Muscles_Emg_Baseline_Integral_Vector(end,:);
    Postural_Muscles_Emg_APA_Integral(i,:) = Postural_Muscles_Emg_APA_Integral_Vector(end,:);
    Postural_Muscles_Emg_CPA_Integral(i,:) = Postural_Muscles_Emg_CPA_Integral_Vector(end,:);


    % Loop to calculate RMS of Baseline, APA, CPA
    for j=1:nc_Postural
        Postural_Muscles_Emg_Baseline_Window_Squared = Postural_Muscles_Emg_Baseline_Window(:,j).^2;
        Postural_Muscles_Emg_APA_Window_Squared = Postural_Muscles_Emg_APA_Window(:,j).^2;
        Postural_Muscles_Emg_CPA_Window_Squared = Postural_Muscles_Emg_CPA_Window(:,j).^2;
        
        for k = 1:(length(Postural_Muscles_Emg_Baseline_Window_Squared)/windowlength_postural_in_lines)
        Postural_Muscles_Emg_Baseline_RMS_Vector(:,j) = (sqrt(mean(Postural_Muscles_Emg_Baseline_Window_Squared(1+windowlength_postural_in_lines*(k-1):windowlength_postural_in_lines*(k)-1))))';
        end
        
        for k = 1:(length(Postural_Muscles_Emg_APA_Window_Squared)/windowlength_postural_in_lines)
        Postural_Muscles_Emg_APA_RMS_Vector(:,j) = (sqrt(mean(Postural_Muscles_Emg_APA_Window_Squared(1+windowlength_postural_in_lines*(k-1):windowlength_postural_in_lines*(k)-1))))';
        Postural_Muscles_Emg_CPA_RMS_Vector(:,j) = (sqrt(mean(Postural_Muscles_Emg_CPA_Window_Squared(1+windowlength_postural_in_lines*(k-1):windowlength_postural_in_lines*(k)-1))))';
        end
     end
    
    Postural_Muscles_Emg_Baseline_RMS(i,:) = Postural_Muscles_Emg_Baseline_RMS_Vector;
    Postural_Muscles_Emg_APA_RMS(i,:) = Postural_Muscles_Emg_APA_RMS_Vector;
    Postural_Muscles_Emg_CPA_RMS(i,:) = Postural_Muscles_Emg_CPA_RMS_Vector;
    
    %Postural_Muscles_Emg_Baseline_RMS_Vector = 0;   %Clearing variables for next iteration
    %Postural_Muscles_Emg_APA_RMS_Vector = 0;
    %Postural_Muscles_Emg_CPA_RMS_Vector = 0;
end

Postural_Muscles_Emg_Integral = [Postural_Muscles_Emg_Baseline_Integral Postural_Muscles_Emg_APA_Integral Postural_Muscles_Emg_CPA_Integral];
Postural_Muscles_Emg_RMS = [Postural_Muscles_Emg_Baseline_RMS Postural_Muscles_Emg_APA_RMS Postural_Muscles_Emg_CPA_RMS];

%%
%%%% FFT ANALYSIS (FREQUENCY DOMAIN)SECTION %%%

fft_signal = fft(All_Channels_Emg_Butter(:,3));
n = length(All_Channels_Emg(:,1));                         % number of samples
f = (0:1/n:1-1/n)*fs;                                      % frequency range
p2 = abs(fft_signal/n);                                    % two-sided power of the DFT
power1 = p2(1:n/2+1);                                       % one-sided power of DFT
power1(2:end-1) = 2*power1(2:end-1);
f1 = fs*(0:(n/2))/n;

%%
%%%% EXPORT DATA SECTION %%%%

Output_File_Path_Full_Integral = [Output_File_Path '\' File_Name '_INTEGRAL' '.xls'];       
Output_File_Path_Full_RMS = [Output_File_Path '\' File_Name '_RMS' '.xls'];       

Repetitions_Row_Label = {'Slow 1º','Slow 2º', 'Slow 3º', 'Fast 1º','Fast 2º', 'Fast 3º', 'Load 1º', 'Load 2º', 'Load 3º'}';
Repetitions_Column_Label = {'Repetition'};
Output_Header_Baseline = repmat({'Baseline'},1,10);
Output_Header_APA = repmat({'APA'},1,10);
Output_Header_CPA = repmat({'CPA'},1,10);

Output_Header_First = [Repetitions_Column_Label Output_Header_Baseline Output_Header_APA Output_Header_CPA];
Output_Header_Second = [Repetitions_Column_Label Postural_Muscles_Name Postural_Muscles_Name Postural_Muscles_Name];
Output_Header = [Output_Header_First; Output_Header_Second];

Output_Data_Integral =  [Repetitions_Row_Label num2cell(Postural_Muscles_Emg_Integral)];
Output_Data_Integral  = [Output_Header; Output_Data_Integral];

Output_Data_RMS =  [Repetitions_Row_Label num2cell(Postural_Muscles_Emg_RMS)];
Output_Data_RMS  = [Output_Header; Output_Data_RMS];

xlswrite(Output_File_Path_Full_Integral, Output_Data_Integral);      
xlswrite(Output_File_Path_Full_RMS, Output_Data_RMS);


%%
%%%% GRAPHICS SECTION %%%%


figure (3)
title ('EMG Activity RIGHT SIDE - PassBand Butter');
xlabel('Time (s)');
ylabel ('Voltage (mV)');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,5),'g');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,11),'b');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,7),'k');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,9),'r');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,13),'c');
hold on
plot (Time_After_Sync(:,:),((All_Channels_Emg_Butter(:,3)/100)+500),'m');
hold on
legend('RAD', 'MULTID','OID','ILEOD','LONGD','DMD')


figure (4)
title ('EMG Activity LEFT SIDE - PassBand Butter');
xlabel('Time (s)');
ylabel ('Voltage (mV)');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,6),'g');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,12),'b');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,15),'k');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,10),'r');
hold on
plot (Time_After_Sync(:,:), All_Channels_Emg_Butter(:,14),'c');
hold on
plot (Time_After_Sync(:,:), ((All_Channels_Emg_Butter(:,3)/100)+500),'m');
hold on
legend('RAE', 'MULTIE','OIE','ILEOE','LONGE','DMD')


figure (5)
title ('EMG Activity RIGHT SIDE - LowPass Butter');
xlabel('Time (s)');
ylabel ('Voltage (mV)');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,1),'g');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,6),'b');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,3),'k');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,4),'r');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,8),'c'); 
hold on
plot (Time_After_Sync(:,:), ((All_Channels_Emg_Butter(:,3)/100)+500),'m');
hold on
legend('RAD', 'MULTID','OID','ILEOD','LONGD','DMD')


figure (6)
title ('EMG Activity LEFT SIDE - LowPass Butter');
xlabel('Time (s)');
ylabel ('Voltage (mV)');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,2),'g');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,7),'b');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,10),'k');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,5),'r');
hold on
plot (Time_After_Sync(:,:), Postural_Muscles_Emg_Butter(:,9),'c'); 
hold on
plot (Time_After_Sync(:,:), ((All_Channels_Emg_Butter(:,3)/100)+500),'m');
hold on
legend('RAE', 'MULTIE','OIE','ILEOE','LONGE','DMD')


figure(7)
plot(f1,power1)
title ('FFT Power Spectrum');
xlabel('Frequency')
ylabel('Power')
