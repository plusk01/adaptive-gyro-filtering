%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Notch Filtering via FFT
%
% Uses a windowed FFT to find the largest frequency component in a
% user-defined range. Peak detection is implemented using a weighted
% average and is smoothed using an alpha filter. A biquad (2nd order IIR)
% notch filter is updated with the detected center frequency. In fact, two
% biquad notch filters are designed with a user-defined separation to
% better notch out the motor noise while keeping phase lag small.
%
% Inspired by Betaflight gyroanalyse.c
%
% CSV data captured in flight from Snapdragon Flight Pro using the
% sensor_imu_tester command, with IMU sampled at 500 Hz.
% 
% Parker Lusk
% 6 Dec 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc;
%% Load raw IMU Data
f = '../data/rawimu-hover.txt';
f = '../data/rawimu-motion.txt';
D = csvread(f,2);

% unpack for convenience
% sequence_number, timestamp(us), updated_timestamp(us),
%   acc_x, acc_y, acc_z, ang_x, ang_y, ang_z, temperature_C 
t = D(:,2);
acc = D(:,4:6);
gyr = D(:,7:9);

% which axis will we inspect?
axis = 1;

% sampling parameters
Fs = 1e6/mean(diff(t));
Ts = 1/Fs;
N = length(t);
n = 0:N-1;
t = Ts*n;

%% Initial LPF'ing

lpfFc = 90;
[lpf.b,lpf.a] = butter(1,lpfFc/(Fs/2));
gyrlpfd = filter(lpf.b,lpf.a,gyr(:,axis));

%% Parameter Setup

% Perform filtering at gyro looprate (i.e., pidDenom = 1)
filterFs = Fs;

FFT_WINDOW_SIZE = 128;
FFT_BIN_COUNT = FFT_WINDOW_SIZE / 2;
DYN_NOTCH_SMOOTH_HZ = 4;
DYN_NOTCH_CALC_TICKS = 3 * 4;

dyn_notch_width_percent = 2;
dyn_notch_q = 360;
dyn_notch_min_hz = 60;
dyn_notch_max_hz = 200;

%% gyroDataAnalyseInit

dynNotch1Ctr = 1 - dyn_notch_width_percent / 100;
dynNotch2Ctr = 1 + dyn_notch_width_percent / 100;
dynNotchQ = dyn_notch_q / 100;
dynNotchMinHz = dyn_notch_min_hz;
dynNotchMaxHz = max(2 * dynNotchMinHz, dyn_notch_max_hz);

if dyn_notch_width_percent == 0
    dualNotch = false;
else
    dualNotch = true;
end

% Notice how fftSamplingRateHz >= 2 * dynNotchMaxHz (i.e., Nyquist if
% dynNotchMaxHz is the highest freq we care about).
gyroLoopRateHz = round(filterFs);
samples = fix(max(1, gyroLoopRateHz / (2 * dynNotchMaxHz)));
fftSamplingRateHz = fix(gyroLoopRateHz / samples);

fftResolution = fftSamplingRateHz / FFT_WINDOW_SIZE;
fftStartBin = fix(max(2, dynNotchMinHz / round(fftResolution)));
% if fftStartBin == Inf, fftStartBin = 2; end
smoothFactor = 2 * pi * DYN_NOTCH_SMOOTH_HZ / (gyroLoopRateHz / 12);

hannWindow = zeros(FFT_WINDOW_SIZE,1);
for i = 0:(FFT_WINDOW_SIZE-1)
    hannWindow(i+1) = (0.5 - 0.5 * cos(2 * pi * i / (FFT_WINDOW_SIZE - 1)));
end

%% gyroDataAnalyseStateInit
maxSampleCount = samples;
state.centerFreq = repmat(dynNotchMaxHz,3,1); % any init value

%% Initialize other state vars

state.sampleCount = 0;
state.oversampledGyroAccumulator = zeros(3,1);
state.downsampledGyroData = zeros(3, FFT_WINDOW_SIZE);

state.circularBufferIdx = 0;
step.updateTicks = 0;
state.updateStep = 0;

state.fftData = zeros(3,FFT_WINDOW_SIZE);
state.rfftData = zeros(3,FFT_WINDOW_SIZE);

%% Initialize filters

notch1 = biquadNotchInit(state.centerFreq(axis) * dynNotch1Ctr, 1/filterFs, dynNotchQ);
notch2 = biquadNotchInit(state.centerFreq(axis) * dynNotch2Ctr, 1/filterFs, dynNotchQ);

%% Plotting setup

% calculate phyiscal frequencies
% f = (0:FFT_BIN_COUNT)*fftResolution;
f = fftSamplingRateHz/2 * linspace(0, 1, FFT_BIN_COUNT+1);

figure(1), clf;
subplot(311); hold on; grid on; title('Gyro');
hP1 = plot(1:maxSampleCount,1:maxSampleCount);
ylabel('rad/s'); xlabel('Sample idx in ring buffer');
subplot(312); hold on; grid on; title('Spectrum');
hP2 = plot(1:maxSampleCount,1:maxSampleCount,'-*');
hP2max = scatter([],[],'*');
hP2cntr = xline(100,'LineWidth',2);
xlabel('Hz');
subplot(313); hold on; grid on; title('After filtering');
hP3 = plot(1:maxSampleCount,1:maxSampleCount,'-*');
xlabel('Hz');

figure(2), clf; hold on;
h2P = plot(t,gyr(:,axis));
h2Pfilt = plot(0,0);
h2Pline = xline(0,'LineWidth',2);
xlim([0 t(end)]);
xlabel('Time (s)'); ylabel('Gyro (rad/s)');
title(['Raw vs Filtered Gyro (axis ' num2str(axis) ')']);

[Hlpf,~] = freqz(lpf.b,lpf.a,1024,fftSamplingRateHz);
[H1,F] = freqz(notch1.b, notch1.a,1024,fftSamplingRateHz);
[H2,~] = freqz(notch2.b, notch2.a,1024,fftSamplingRateHz);
H = Hlpf .* H1 .* H2;
figure(3), clf;
subplot(311); hold on; grid on; box on; title('Filter Design');
h3P1 = plot(F, 20*log10(abs(H)));
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
subplot(312); hold on; grid on; box on;
h3P2 = plot(F, rad2deg(angle(H)));
xlabel('Frequency (Hz)'); ylabel('Phase (degrees)');
subplot(313); hold on; grid on; box on;
h3P3 = plot(F, -unwrap(angle(H)) ./ (2*pi*F/fftSamplingRateHz));
xlabel('Frequency (Hz)'); ylabel('Sample lag');

gyrofiltered = zeros(1,N);

%% Main loop
for i = 1:N
    k = n(i);
    
    %
    % gyro.c / gyro_filter_impl.c : filterGyro()
    %
    
    % gyroDataAnalysePush
    state.oversampledGyroAccumulator = state.oversampledGyroAccumulator + gyr(i,:)';
    
    % perform notch filtering
    y = gyrlpfd(i);
    [notch1, y] = biquadApplyDF1(notch1, y);
    [notch2, y] = biquadApplyDF1(notch2, y);
    gyrofiltered(i) = y;
    
    %
    % gyro.c / gyroanalyse.c : gyroDataAnalyse()
    %
    
    state.sampleCount = state.sampleCount + 1;
    
    % Downsample to FFT rate via averaging
    if state.sampleCount == maxSampleCount
        state.sampleCount = 0;
        sample = state.oversampledGyroAccumulator / maxSampleCount;
        state.downsampledGyroData(:,state.circularBufferIdx+1) = sample;
        state.oversampledGyroAccumulator = zeros(3,1);
       
        state.circularBufferIdx = mod(state.circularBufferIdx + 1, FFT_WINDOW_SIZE);
        
%         state.updateTicks = DYN_NOTCH_CALC_TICKS;
        state.updateTicks = 1;
    end
    
    % CPU util management by breaking gyro analysis into smaller steps
    if state.updateTicks == 0, continue; end
    state.updateTicks = 0; %state.updateTicks - 1;
    
    %
    % gyroanalyse.c : gyroDataAnalyseUpdate()
    %
    
    updateFFTPlots = false;
    
    STEP_ARM_CFFT_F32 = 0;
    STEP_BITREVERSAL = 1;
    STEP_STAGE_RFFT_F32 = 2;
    STEP_ARM_CMPLX_MAG_F32 = 3;
    STEP_CALC_FREQUENCIES  = 4;
    STEP_UPDATE_FILTERS = 5;
    STEP_HANNING = 6;
    STEP_COUNT = 7;
    
    if state.updateStep == STEP_ARM_CFFT_F32
        state.fftData(1,:) = fft(state.fftData(1,:)) / N;
        state.fftData(2,:) = fft(state.fftData(2,:)) / N;
        state.fftData(3,:) = fft(state.fftData(3,:)) / N;
        
    elseif state.updateStep == STEP_BITREVERSAL
        state.updateStep = state.updateStep + 1;
%     elseif state.updateStep == STEP_STAGE_RFFT_F32
%         state.rfftData(1,:) = fftshift(state.fftData(1,:));
%         state.rfftData(2,:) = fftshift(state.fftData(2,:));
%         state.rfftData(3,:) = fftshift(state.fftData(3,:));
        state.rfftData = state.fftData;
        
    elseif state.updateStep == STEP_ARM_CMPLX_MAG_F32
        state.fftData = abs(state.rfftData);
%         state.fftData = state.fftData(:,FFT_WINDOW_SIZE/2:end);
        state.updateStep = state.updateStep + 1;
%     elseif state.updateStep == STEP_CALC_FREQUENCIES

    dataMax = 0;
    dataMin = 1;
    binMax = 0;
    dataMinHi = 1;
    for ii = fftStartBin:FFT_BIN_COUNT
        if state.fftData(axis,ii+1) > state.fftData(axis,ii) % bin height increased
            if state.fftData(axis,ii+1) > dataMax
                dataMax = state.fftData(axis,ii+1);
                binMax = ii;
            end
        end
    end
    if binMax == 0 % no bin increase, hold prev max bin
        binMax = fix(state.centerFreq(axis) / fftResolution);
    else % there was a max, find min
        for ii = binMax-1:-1:1 % look for min below max
            dataMin = state.fftData(axis,ii+1);
            if state.fftData(axis,ii) > state.fftData(axis,ii+1), break; end
        end
        for ii = binMax+1:(FFT_BIN_COUNT-1) % look for min above max
            dataMinHi = state.fftData(axis,ii+1);
            if state.fftData(axis,ii+1) < state.fftData(axis,ii), break; end
        end
    end
    dataMin = min(dataMin, dataMinHi);
    
    % accumulate fftSum and fftWeightedSum from peak bin, and shoulder bins either side of peak
    squaredData = state.fftData(axis,binMax+1) ^ 2;
    fftSum = squaredData;
    fftWeightedSum = squaredData * binMax;
    
    % accumulate upper shoulder unless it would be FFT_BIN_COUNT
    shoulderBin = binMax + 1;
    if shoulderBin < FFT_BIN_COUNT
        squaredData = state.fftData(axis,shoulderBin+1) ^ 2;
        fftSum = fftSum + squaredData;
        fftWeightedSum = fftWeightedSum + squaredData * shoulderBin;
    end
    
    % accumulate lower shoulder unless lower shoulder would be bin 0 (DC)
    if binMax > 1
        shoulderBin = binMax - 1;
        squaredData = state.fftData(axis,shoulderBin+1) ^ 2;
        fftSum = fftSum + squaredData;
        fftWeightedSum = fftWeightedSum + squaredData * shoulderBin;
    end
    
    % get centerFreq in hz from weighted bins (weighted mean)
    centerFreq = 0;
    fftMeanIndex = 0;
    if fftSum > 0
        fftMeanIndex = fftWeightedSum / fftSum;
        centerFreq = fftMeanIndex * fftResolution;
    else
        centerFreq = state.centerFreq(axis);
    end
    centerFreq = constrain(centerFreq, dynNotchMinHz, dynNotchMaxHz);
    
    % LPF-style dynamic smoothing
%     dynamicFactor = constrain(dataMax / dataMin, 1, 2);
    dynamicFactor = 1;
    smoothFactor = 0.1;
    state.centerFreq(axis) = state.centerFreq(axis) + smoothFactor * dynamicFactor * (centerFreq - state.centerFreq(axis));

    updateFFTPlots = true;

    elseif state.updateStep == STEP_UPDATE_FILTERS
        notch1 = biquadNotchUpdate(notch1, state.centerFreq(axis) * dynNotch1Ctr, 1/filterFs, dynNotchQ);
        notch2 = biquadNotchUpdate(notch2, state.centerFreq(axis) * dynNotch2Ctr, 1/filterFs, dynNotchQ);
        state.updateStep = state.updateStep + 1;
%     elseif state.updateStep == STEP_HANNING
        ringBufIdx = FFT_WINDOW_SIZE - state.circularBufferIdx;
        state.fftData(:,1:ringBufIdx) = state.downsampledGyroData(:,(state.circularBufferIdx+1):FFT_WINDOW_SIZE) .* repmat(hannWindow(1:ringBufIdx)',3,1);
        if state.circularBufferIdx > 0
            state.fftData(:,(ringBufIdx+1):FFT_WINDOW_SIZE) = state.downsampledGyroData(:,1:state.circularBufferIdx) .* repmat(hannWindow(ringBufIdx+1:FFT_WINDOW_SIZE)',3,1);
        end
        
        % make copy of buffer data just for extra analysis
        state.data = state.fftData;
    end
    
    state.updateStep = mod(state.updateStep + 1, STEP_COUNT);
    
    %
    % Visualization
    %

    if updateFFTPlots
        set(hP1,'XData',1:FFT_WINDOW_SIZE,'YData',state.downsampledGyroData(axis,:));
        set(hP2,'XData',f,'YData',state.fftData(axis,1:FFT_BIN_COUNT+1));
        set(hP2max,'XData',f(binMax+1),'YData',state.fftData(axis,binMax+1));
        set(hP2cntr,'Value',state.centerFreq(axis));
        figure(1); subplot(312); title(['centerFreq = ' num2str(round(state.centerFreq(axis))) ' Hz'])

        s = constrain(i-FFT_WINDOW_SIZE, 1, i);
        Y = fft(gyrofiltered(s:i),FFT_WINDOW_SIZE)/N;
        set(hP3,'XData',f,'YData',abs(Y(1:FFT_BIN_COUNT+1)));
        drawnow
    end
    
    set(h2Pfilt,'XData',t(1:i),'YData',gyrofiltered(1:i));
    set(h2Pline,'Value',t(i));
    
    [H1,F] = freqz(notch1.b, notch1.a,1024,fftSamplingRateHz);
    [H2,~] = freqz(notch2.b, notch2.a,1024,fftSamplingRateHz);
    H = Hlpf .* H1 .* H2;
    set(h3P1,'YData',20*log10(abs(H)));
    set(h3P2,'YData',rad2deg(angle(H)));
    set(h3P3,'YData',-unwrap(angle(H)) ./ (2*pi*F/fftSamplingRateHz));
end

%% Helpers
function x = constrain(x, low, high)
    if x < low, x = low; end
    if x > high, x = high; end
end

function filter = biquadNotchInit(fc, Ts, Q)

    % normalized frequency in [0, 2pi]
    omega = 2 * pi * fc * Ts;
    sn = sin(omega);
    cs = cos(omega);
    alpha = sn / (2 * Q);
    
    filter.fc = fc;
    filter.Ts = Ts;
    filter.Q = Q;
    
    % notch (b num, a denom)
    b0  = 1;
    b1 = -2 * cs;
    b2 = 1;
    a0 = 1 + alpha;
    a1 = -2 * cs;
    a2 = 1 - alpha;
    
    % normalize into biquad form (a0 == 1)
    filter.b0 = b0 / a0;
    filter.b1 = b1 / a0;
    filter.b2 = b2 / a0;
    filter.a1 = a1 / a0;
    filter.a2 = a2 / a0;
    
    % convenience:
    filter.b = [filter.b0 filter.b1 filter.b2];
    filter.a = [1 filter.a1 filter.a2];
    
    % initialize delay elements
    filter.x1 = 0;
    filter.x2 = 0;
    filter.y1 = 0;
    filter.y2 = 0;
end

% requires a direct form 1 (DF1) apply implementation to allow changing coeffs
function filter = biquadNotchUpdate(filter, fc, Ts, Q)
    % backup state
    x1 = filter.x1;
    x2 = filter.x2;
    y1 = filter.y1;
    y2 = filter.y2;
    
    filter = biquadNotchInit(fc, Ts, Q);
    
    % restore state
    filter.x1 = x1;
    filter.x2 = x2;
    filter.y1 = y1;
    filter.y2 = y2;
end

% slightly less precise than DF2, but works in dynamic mode
function [filter, y] = biquadApplyDF1(filter, x)
    y = filter.b0 * x + filter.b1 * filter.x1 + filter.b2 * filter.x2 - (filter.a1 * filter.y1 + filter.a2 * filter.y2);
    
    % shift feedback delay lines
    filter.x2 = filter.x1;
    filter.x1 = x;
    
    % shift feedforward delay lines
    filter.y2 = filter.y1;
    filter.y1 = y;
end