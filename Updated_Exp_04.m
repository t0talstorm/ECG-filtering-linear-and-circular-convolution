%% ECG FILTERING: LINEAR AND CIRCULAR CONVOLUTION
% MATLAB R2026a 
%
% Healthcare application:
% Comparison of linear and circular convolution in ECG filtering.
%
% x[n] = Noisy ECG signal
% h[n] = Five-sample moving-average FIR filter
%
% Linear convolution preserves the actual ECG time sequence.
% Circular convolution may wrap the ending samples to the beginning.

clc;
clear;
close all;

%% 1. CREATE A SIMPLE ECG SIGNAL

Fs = 200;                     % Sampling frequency
duration = 3;                 % ECG duration in seconds

t = 0:1/Fs:duration-1/Fs;
N = length(t);

% Heartbeat locations
% Last heartbeat is placed close to the end
beatTimes = [0.5 1.3 2.1 2.97];

ecgClean = zeros(size(t));

for beat = beatTimes

    % P wave
    P = 0.12*exp(-((t-(beat-0.18)).^2)/(2*0.025^2));

    % Q wave
    Q = -0.15*exp(-((t-(beat-0.03)).^2)/(2*0.010^2));

    % R wave
    R = 1.00*exp(-((t-beat).^2)/(2*0.012^2));

    % S wave
    S = -0.25*exp(-((t-(beat+0.03)).^2)/(2*0.012^2));

    % T wave
    TWave = 0.25*exp(-((t-(beat+0.22)).^2)/(2*0.050^2));

    ecgClean = ecgClean + P + Q + R + S + TWave;
end

%% 2. ADD 50 Hz POWER-LINE INTERFERENCE

noise = 0.08*sin(2*pi*50*t);

% x[n]: Noisy ECG signal
x = ecgClean + noise;

%% 3. DEFINE FIVE-SAMPLE FIR FILTER

M = 5;

h = ones(1,M)/M;

% h[n] = [0.2 0.2 0.2 0.2 0.2]

%% 4. LINEAR CONVOLUTION

yLinearFull = conv(x,h,'full');

tLinear = (0:length(yLinearFull)-1)/Fs;

% First N samples for comparison
yLinear = yLinearFull(1:N);

%% 5. CIRCULAR CONVOLUTION

yCircular = real(ifft(fft(x,N).*fft(h,N)));

%% 6. ZERO-PADDED FFT CONVOLUTION

Nfft = length(x)+length(h)-1;

yFFTLinear = real(ifft( ...
    fft(x,Nfft).*fft(h,Nfft)));

%% 7. DISPLAY RESULTS

verificationError = max(abs(yLinearFull-yFFTLinear));

wrapError = yCircular-yLinear;

fprintf('\nECG CONVOLUTION ANALYSIS\n');
fprintf('------------------------------------------\n');

fprintf('Number of ECG samples      = %d\n',N);
fprintf('Filter length              = %d\n',M);
fprintf('Linear output length       = %d\n',length(yLinearFull));
fprintf('Circular output length     = %d\n',length(yCircular));

fprintf('Zero-padded FFT error      = %.12f\n', ...
        verificationError);

fprintf('Maximum wrap-around error  = %.6f\n', ...
        max(abs(wrapError)));

%% 8. PLOT THE RESULTS

figure('Color','white', ...
       'Name','ECG Convolution Comparison', ...
       'NumberTitle','off', ...
       'Position',[80 40 1400 850]);

tiledlayout(4,1, ...
            'TileSpacing','compact', ...
            'Padding','compact');

%% Plot 1: Noisy and clean ECG

nexttile;

p1 = plot(t,x,'r','LineWidth',1.2);
hold on;

p2 = plot(t,ecgClean,'b','LineWidth',1.6);

grid on;
box on;

set(gca, ...
    'Color','white', ...
    'XColor','black', ...
    'YColor','black', ...
    'GridColor',[0.75 0.75 0.75], ...
    'FontWeight','bold');

xlabel('Time (seconds)','FontWeight','bold');
ylabel('ECG Amplitude','FontWeight','bold');

title('x[n]: ECG Signal with 50 Hz Interference', ...
      'Color','black', ...
      'FontWeight','bold');

lgd1 = legend([p1 p2], ...
              'Noisy ECG', ...
              'Clean ECG', ...
              'Location','northeast');

lgd1.Color = 'white';
lgd1.TextColor = 'black';
lgd1.EdgeColor = 'black';
lgd1.FontWeight = 'bold';

%% Plot 2: FIR filter h[n]

nexttile;

stem(0:M-1,h,'g','filled','LineWidth',1.5);

grid on;
box on;

set(gca, ...
    'Color','white', ...
    'XColor','black', ...
    'YColor','black', ...
    'GridColor',[0.75 0.75 0.75], ...
    'FontWeight','bold');

xlabel('Sample Index, n','FontWeight','bold');
ylabel('Coefficient','FontWeight','bold');

title('h[n]: Five-Sample Moving-Average FIR Filter', ...
      'Color','black', ...
      'FontWeight','bold');

ylim([0 0.22]);

%% Plot 3: Linear convolution

nexttile;

plot(tLinear,yLinearFull,'b','LineWidth',1.6);
hold on;

xline(duration,'k--', ...
      'End of ECG Recording', ...
      'LineWidth',1.3, ...
      'FontWeight','bold', ...
      'LabelVerticalAlignment','middle');

grid on;
box on;

set(gca, ...
    'Color','white', ...
    'XColor','black', ...
    'YColor','black', ...
    'GridColor',[0.75 0.75 0.75], ...
    'FontWeight','bold');

xlabel('Time (seconds)','FontWeight','bold');
ylabel('ECG Amplitude','FontWeight','bold');

title('Linear Convolution: Filter Response Continues After the End', ...
      'Color','black', ...
      'FontWeight','bold');

%% Plot 4: Linear and circular comparison

nexttile;

p3 = plot(t,yLinear,'b','LineWidth',2);
hold on;

p4 = plot(t,yCircular,'r--','LineWidth',1.6);

grid on;
box on;

set(gca, ...
    'Color','white', ...
    'XColor','black', ...
    'YColor','black', ...
    'GridColor',[0.75 0.75 0.75], ...
    'FontWeight','bold');

xlabel('Time (seconds)','FontWeight','bold');
ylabel('ECG Amplitude','FontWeight','bold');

title('Circular Convolution Produces Wrap-Around at the Beginning', ...
      'Color','black', ...
      'FontWeight','bold');

lgd2 = legend([p3 p4], ...
              'Linear Convolution', ...
              'Circular Convolution', ...
              'Location','northwest');

lgd2.Color = 'white';
lgd2.TextColor = 'black';
lgd2.EdgeColor = 'black';
lgd2.FontWeight = 'bold';

sgtitle('Why Linear Convolution Is Preferred for ECG Filtering', ...
        'Color','black', ...
        'FontWeight','bold', ...
        'FontSize',16);

%% 9. ZOOMED VIEW OF THE BEGINNING

figure('Color','white', ...
       'Name','ECG Boundary Artefact', ...
       'NumberTitle','off', ...
       'Position',[150 120 1000 520]);

p5 = plot(t,yLinear,'b','LineWidth',2);
hold on;

p6 = plot(t,yCircular,'r--','LineWidth',1.8);

grid on;
box on;

set(gca, ...
    'Color','white', ...
    'XColor','black', ...
    'YColor','black', ...
    'GridColor',[0.75 0.75 0.75], ...
    'FontWeight','bold');

xlim([0 0.15]);

xlabel('Time (seconds)','FontWeight','bold');
ylabel('ECG Amplitude','FontWeight','bold');

title('Artificial Beginning Created by Circular Wrap-Around', ...
      'Color','black', ...
      'FontWeight','bold');

lgd3 = legend([p5 p6], ...
              'Linear: Correct Beginning', ...
              'Circular: Wrapped Final Samples', ...
              'Location','northeast');

lgd3.Color = 'white';
lgd3.TextColor = 'black';
lgd3.EdgeColor = 'black';
lgd3.FontWeight = 'bold';