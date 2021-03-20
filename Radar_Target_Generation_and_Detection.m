clear all
clc;

%% Radar Specifications 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency of operation = 77 GHz  %
% Max Range = 200 m                %
% Range Resolution = 1 m           %
% Max Velocity = 100 m/s           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maxRange = 200;
rangeResolution = 1;
maxV = 100;
c = 3e8;

% User defined initial position (m) and velocity (m/s) of target. Note : Velocity remains contant
R = 110;
V = -20;


%% FMCW Waveform Generation

% Bandwidth
B = c / (2 * rangeResolution);

% Chirp Time
Tchirp = 5.5 * (2 * maxRange / c);

% Slope of the FMCW chirp
slope = B / Tchirp;

%Operating carrier frequency of Radar 
fc= 77e9;
                                                       
%The number of chirps in one sequence.
Nd = 128;

%The number of samples on each chirp. 
Nr = 1024;

% Timestamp for running the displacement scenario for every sample on each chirp
t = linspace(0, Nd*Tchirp, Nr*Nd);


%% Signal generation and Moving Target simulation
% Running the radar scenario over the time. 

% Range covered
r_t = R + V * t;

% Time delay
td = 2 * r_t / c;

% Transmitted signal
Tx = cos(2*pi*(fc*t+slope*t.^2/2));

% Received signal
Rx = cos(2*pi*(fc*(t-td)+slope*(t-td).^2/2));

% Beat signal
Mix = Tx .* Rx;


%% RANGE MEASUREMENT

Mix = reshape(Mix, [Nr,Nd]);
sig_fft = fft(Mix, Nr) ./ Nr;
sig_abs = abs(sig_fft);
sig_fft_clip = sig_abs(1:Nr/2);

figure(1)
plot(sig_fft_clip)
title('Range from First FFT')
xlabel('Range (m)')
axis ([0 200 0 1]);



%% RANGE DOPPLER RESPONSE

Mix = reshape(Mix, [Nr,Nd]);
sig_fft2 = fft2(Mix,Nr,Nd);
sig_fft2 = sig_fft2(1:Nr/2,1:Nd);
sig_fft2 = fftshift (sig_fft2);
RDM = abs(sig_fft2);
RDM = 10 * log10(RDM) ;

doppler_axis = linspace(-100,100,Nd);
range_axis = linspace(-200,200,Nr/2)*((Nr/2)/400);
figure(2)
surf(doppler_axis,range_axis,RDM);
title('Range-Doppler Map')
xlabel('Doppler Velocity (m/s)')
ylabel('Range (m)')


%% CFAR Implementation

% Number of Training Cells in both the dimensions.
Tr = 12;
Td = 10;

% Number of Guard Cells in both dimensions.
Gr = 5;
Gd = 5;

% Offset the threshold by SNR value in dB
offset = 6;


newRDM = zeros(size(RDM));
for i = Tr + Gr + 1 : Nr / 2 - (Tr + Gr)
    for j = Td + Gd + 1 : Nd - (Td + Gd)
        noise_level = zeros(1,1);
        for p = i - (Tr + Gr) : i + (Tr + Gr)
            for q = j - (Td + Gd) : j + (Td + Gd)
                
                if (abs(i - p) > Gr || abs(j - q) > Gd)
                    noise_level = noise_level + db2pow(RDM(p,q));
                end
                
            end
            
        end
        
        cellSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1)-(2*Gr+1)*(2*Gd+1);
        noise_level = pow2db(noise_level / cellSize);
        threshold = noise_level + offset;
        
        CUT = RDM(i,j);
        if CUT > threshold
            newRDM(i,j) = 1;
        end
        
    end
    
end


figure(3)
surf(doppler_axis,range_axis, newRDM);
colorbar;
title('CFAR Implementation');
xlabel('Doppler Velocity (m/s)')
ylabel('Range (m)')
