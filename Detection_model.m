%% SECTION TITLE
% DESCRIPTIVE DTMF decoder
close all
clear all
clc

% Read the audio file
[y,Fs] = audioread('0123456789.wav');
% [y,Fs] = audioread('numero2.wav');
% [y,Fs] = audioread('numero3.wav');
% [y,Fs] = audioread('numero4.wav');

% Ts is the sampling period
Ts = 1/Fs;

npts = length(y);
% time vector
t = (0:1:npts-1)*Ts;
fig=1;

%% Plot Figure

figure(fig)
plot(t,y)
xlabel('Time in seconds')
ylabel('Sampled Siignals')
title('Plot of Original Signal in TIme Domain')

%% Spectrum Analysis of the signal
%M = length(y);
M = 256;

f1=linspace(0,Fs,M+1);
f1=f1(1:end-1);
f2=f1/Fs;
f1=f1-Fs/2;

Y = fft(y,M);

fig=fig+1;
figure(fig)
plot(f2,abs(Y))
title('Spectrum of the signal before filtering')
xlabel('Frequency in Hz')
ylabel('Spectrum of the signal')

%% Spectrum analysis of each window in 65 ms
D = 65e-3;
nn = round(D/Ts);

% Reshaping the signal in intervals of the window size
n = length(y);
fint = 1/n;
f = 0:fint:(n-1)*fint;

nb = ceil(length(y)/nn)*nn - length(y);
ye1 = [y.' zeros(1,nb)];
y_tab = reshape(ye1', nn, [])';
y_tab2 = zeros(size(y_tab));

t_tab1 = [t  zeros(1,nb)];
t_tab2 = reshape(t_tab1', nn, [])';

%%

% minimum power initialization
power_min = 20;


%% Frequency vector
n = length(y);
fint = 1/n;
f = 0:fint:(n-1)*fint;

% Initialization of the vector
Tone = zeros(size(y_tab,1),1);

for k = 1:size(y_tab,1)
    
    power = sum(y_tab(k,:).^2)/nn;% power of the signal
    power_dB = (10*log10(power/10.^-3)); %power of the signal in DB
    
    Tone(k) = power_dB > power_min; % to select all windows with power greater 20dB
    
end

s1 = sum(Tone);

power = zeros(size(y_tab,1),1);
power_dB = zeros(size(y_tab,1),1);
Tone2 = zeros(s1,2);

%% parameters for cheby1 () filter function
aa = 7;
bb = 10;
cc = 0.55;
dd = 0.93;

for k = 1:size(y_tab,1)
    
    if  k == 1      %
        
        Y_TAB1 = fft(y_tab(k,:),n);
        
        Y_TAB1_mag = abs(Y_TAB1);
        
        [peak1, loc1] = findpeaks(Y_TAB1_mag(1:(n+1)/2),'SortStr','descend');
        Tone2(k,1) = f(loc1(1)); % to get the location of the first frequency
        N1 = Tone2(k,1);
        
        %[m1, n1] = max(Y_TAB1_mag(1:(n+1)/2));
        %Tone2(k,1) = f(n1-1);
        %N1 = Tone2(k,1);
        
        [b,a] = cheby1(aa,bb,[cc*N1 dd*N1],'stop');
        
        y_tab2(k,:) = filter(b,a,y_tab(k,:));
        Y_TAB2=fft(y_tab2(k,:),n);
        Y_TAB2_mag = abs(Y_TAB2);
        
        [peak2, loc2] = findpeaks(Y_TAB2_mag(1:(n+1)/2),'SortStr','descend');
        
        Tone2(k,2) = f(loc2(1));
        
    else
        
        power(k-1) = sum(y_tab(k-1,:).^2)/nn;% power of the signal
        power_dB(k-1) = (10*log10(power(k-1)/10.^-3)); %power of the signal in DB
        
        power(k) = sum(y_tab(k,:).^2)/nn;% power of the signal
        power_dB(k) = (10*log10(power(k)/10.^-3)); %power of the signal in DB
        
        if power_dB(k-1) < power_min && power_dB(k) > power_min % to avoid taking duplicate frequencies
            
            Y_TAB1 = fft(y_tab(k,:),n); %to get the spectrum
            
            Y_TAB1_mag = abs(Y_TAB1); %to get the magniure
            
            [peak1, loc1] = findpeaks(Y_TAB1_mag(1:(n+1)/2),'SortStr','descend');
            
            %[m1, n1] = max(Y_TAB1_mag(1:(n+1)/2));
            %Tone2(k,1) = f(n1);
            Tone2(k,1) = f(loc1(1));
            %N1 = loc1(1)/Fs;
            
            N1 = Tone2(k,1);
            [b,a] = cheby1(aa,bb,[cc*N1 dd*N1],'stop');% to design the filter
            
            y_tab2(k,:) = filter(b,a,y_tab(k,:)); % to apply the filter 
            Y_TAB2=fft(y_tab2(k,:),n); %to get the spectrum
            Y_TAB2_mag = abs(Y_TAB2); % to get the magnitude
            
            [peak2, loc2] = findpeaks(Y_TAB2_mag(1:(n+1)/2),'SortStr','descend');
            Tone2(k,2) = f(loc2(1)); % to get the location of the second frequency
            
            %[m2, n2] = max(Y_TAB1_mag(1:(n+1)/2));
            %Tone2(k,2) = f(n2);
            
        end
    end
end
% the first and second column were sorted to make the first column to be
% the low frequency and second column the high frequency
Tone3 =reshape(nonzeros(Tone2),[],2)*Fs;
Tone4 = sort(Tone3.','ascend');
Tone5 = Tone4';
number = zeros(length(Tone5),1);

%% To compare the frequencies in the signal to the DTMF, considering the 1.5% frequency tolerance

for k = 1:length(Tone5)
    
    if Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.455 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = 1;
        
    elseif  Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.455 && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 2;
        
    elseif Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.455 && Tone5(k,2) >= 1454.845 && Tone5(k,2) <= 1499.155
        number(k) = 3;
        
    elseif  Tone5(k,1) >= 686.545 && Tone5(k,1) <= 707.455 && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'A';
        
        
        
        
    elseif Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = 4;
        
    elseif  Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55  && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 5;
        
    elseif Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55  && Tone5(k,2) >= 1454.845 && Tone5(k,2) <= 1499.155
        number(k)= 6;
        
    elseif  Tone5(k,1) >= 758.45 && Tone5(k,1) <= 781.55  && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'B';
        
        
        
        
    elseif Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = 7;
        
    elseif  Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78  && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 8;
        
    elseif Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78  && Tone5(k,2) >= 1454.845 && Tone5(k,2) <= 1499.155
        number(k) = 9;
        
    elseif  Tone5(k,1) >= 839.22 && Tone5(k,1) <= 864.78  && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'C';
        
        
        
        
    elseif Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115 && Tone5(k,2) >= 1190.865 && Tone5(k,2) <= 1227.135
        number(k) = '*';
        
    elseif  Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115  && Tone5(k,2) >= 1315.96 && Tone5(k,2) <= 1356.04
        number(k) = 0;
        
    elseif Tone5(k,1) >= 926.885 && Tone5(k,1) <= 955.115  && Tone5(k,2) >= 1454.845 && Tone5(k,2) <= 1499.155
        number(k)= '#';
        
    elseif  Tone5(k,1) >= 926.885  && Tone5(k,1) <= 955.115 && Tone5(k,2) >= 1612.445 && Tone5(k,2) <= 1661.555
        number(k) = 'D';
        
    end
end

number'

