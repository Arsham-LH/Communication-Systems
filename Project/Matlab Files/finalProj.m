%% record voice
continuous_freq = 44100;
recObj = audiorecorder(continuous_freq,16,1) ;
recDuration = 5 ; % duration of voice in seconds
disp("Begin Speaking");
recordblocking(recObj , recDuration );
disp("End of Recording");
sig = getaudiodata(recObj);
Fsig = recObj.SampleRate;
%% save voice
save('voice','sig','Fsig');
%% load voice
voice = load("voice.mat");
Fsig = voice.Fsig;
sig = voice.sig;
%% play voice (method 1)
sound(sig,Fsig);
%% play voice (method 2)
ap = audioplayer(sig , Fsig);
play(ap);
%% simultaneous

%% apply the system
% first record a voice and put it in sig ( run the first two sections or the third section)
Fs_system = continuous_freq/10;
v = 8;
A = 1;
continuous_freq = 44100 ;
% always make sure : v*Fs_system < continuous_freq
N0 = 10^(-7);
B = 20000;
r = v*Fs_system;
beta = r/2;
[t,outSignal,bi1,x2,x3,bi2] = main(sig , Fsig , Fs_system , v , r, beta , A , continuous_freq , N0 , B);
Fout = continuous_freq ;
d = bi1-bi2; c = d(d~=0); ber = length(c)/length(bi1);
%% play
sound(outSignal,Fout);
save('outvoice_fs_4410_B_20000_beta_halfr.mat','outSignal','Fout');
%%
A = [1 5 10 15 20];
continuous_freq = 44100 ;
Fs_system = continuous_freq/12;
v = 8;
r = v*Fs_system;
beta = r/2;
B = 5000;
ber = zeros(5,20);
n=(0:5:95)*10^(-7);
for i=4:1:5
    for j=1:1:20
        [t,outSignal,bi1,x2,x3,bi2] = main(sig , Fsig , Fs_system , v , r, beta , A(i) , continuous_freq , n(j) , B);
        Fout = continuous_freq ;
        d = bi1-bi2; c = d(d~=0); ber(i,j) = length(c)/length(bi1);
    end
end



%% Plotting figures
figure;
subplot(1,2,1)
plot(t,sig);
title('Input');
xlabel('t(sec)');
subplot(1,2,2)
plot(t,outSignal);


%% 
sigma2 = n*B ;
figure;
p1 = plot(sigma2(1:20) , ber(1,:)); hold on
p2 = plot(sigma2(1:20) , ber(2,:));
p3 = plot(sigma2(1:20) , ber(3,:));
p4 = plot(sigma2(1:20) , ber(4,:));
p5 = plot(sigma2(1:20) , ber(5,:));
legend(["A="+A(1),"A="+A(2),"A="+A(3),"A="+A(4),"A="+A(5)]);
%% Functions

function [t,outputSignal,x1,x2,x3,x4] = main(inputSignal , inputSampleRate , requiredSampleRate , numberOfBits , r, beta , A , continuous_fs , N0 , B)
%     r = numberOfBits*requiredSampleRate ;
    x1 = adc(inputSignal , inputSampleRate , requiredSampleRate , numberOfBits );
    [t,x2] = linecoder(x1 , beta , r , A , continuous_fs);
    x3 = channel(x2 , N0 , B , continuous_fs);
    x4 = linedecoder(x3 , r , continuous_fs );
    outputSignal = dac(x4 , numberOfBits , requiredSampleRate , continuous_fs);
end

function y = adc( inputSignal , inputSampleRate , requiredSampleRate , numberOfBits )
    x = changeSampleRate(inputSignal , inputSampleRate , requiredSampleRate);
    y = digitalize( x , numberOfBits);
end
function y = changeSampleRate( signal , f1 , f2) % change sample rate from f1 to f2
    L = length(signal);
    T = L/f1 ; % duration of signal
    Ly = round(T*f2) ; % length of y
    t = round((0:1:Ly-1)*f1/f2)+1; % scale index
    y = signal(t);
    

end
function y = digitalize( signal , numberOfBits) % find binary code for one quantized number
    x = quantize( signal , numberOfBits);
    z = de2bi( x , numberOfBits );
    y = reshape(z.' , length(signal)*numberOfBits , 1);
end
function y = quantize( signal , numberOfBits )
    n = 2^(numberOfBits-1) ; % n = number of codes / 2
    x = changeAmplitude( signal )*n; % x => (-n , n)
    y = (floor(x)+n); % x => [0 , 1 , ... , 2n-1]
end
function y = changeAmplitude( x ) % change amplitude of x to a little less than 1
    y = x/(max(abs(x),[],'all')*1.0001);
end

function [t,lineCode_x] = linecoder(signal, beta,r,A,continuous_fs)
    %*****better to choose r as a fraction of system_fs, for example r=system_fs/1000******
    D = 1/r;
    L = length(signal);
    T = D*L;
    t = 1/continuous_fs:1/continuous_fs:T; %WARNING: THIS t IS DIFFERENT FROM t OF THE MAIN SIGNAL!!!
    t_len = length(t);
    t_posNeg=cat(2,-fliplr(t),0,t); %t_posNeg contains time from -t to +t
    
    lineCode_x_ampl = (signal-0.5)*2*A; %creating amplitude +-A from bits 1&0
    
     
    lineCode_x_zeroIntrp = zeros(t_len+1,1);
    %lineCode_x_zeroIntrp(1:t_D_len:end) = lineCode_x_ampl; %zero interpolation for convolving
    for i=1:1:L
        lineCode_x_zeroIntrp(round((i-1)*continuous_fs/r)+1)=lineCode_x_ampl(i);
    end
    
    p = cospi(2*beta*t_posNeg) ./ (1-(4*beta*t_posNeg).^2) .* sinc(r*t_posNeg); %p(t), a column vector. length = 2*length(t)+1 = length(t_posNeg)
    p=p.'; %converting to row vector

    %lineCode_x = zeros((encode_vect_len-1)*t_D_len+1,1);
    %disp("line 190");
    lineCode_x = conv(lineCode_x_zeroIntrp,p);
    lineCode_x = lineCode_x(t_len+1:2*t_len); %removing negative parts and last part (that exceeds t_len)

    %******* at the end, sample 1 from encode_vect corresponds to sample 1
    %from lineCode_x (or t), sample 2 corresponds to sample 101 from
    %lineCode_x,... and sample i corresponds to i*(t_D_len-1)+1 **********
    
    %*** lineCode_x starts with the first bit peak and finishes with last bit peak ***
end

function y = channel(signal , N0 , B , fs)
    noiseVar = N0*B ;
    noise = wgn(size(signal,1),size(signal,2),noiseVar,'linear');
    % noise = sqrt(noiseVar)*randn(size(signal,1),size(signal,2));
    y = signal + noise;
    L = length(y);
    if fs>B 
        f = fft(y);
        n = round(B/(fs/2)*L);
        f(n+2:L-n) = 0;
        y = ifft(f);
    end
end

function y = linedecoder( signal , r , continuous_fs )
    x = changeSampleRate(signal , continuous_fs , r); % sampling
    L = length(x);
    y = zeros(L,1);
    y(x>0)=1; % binary signal
end

function z = dac(signal , numberOfBits , Fsignal , continuous_f)
    L = length(signal);
    Ly = round(L/numberOfBits);
    x = reshape(signal , numberOfBits , []); % divide symbols
    y = bi2de(x.')+(-2^(numberOfBits-1)+0.5); % decode binary to quantized numbers
    z = increaseSampleRate( y , Fsignal , continuous_f);
    z = z/max(abs(z),[],'all');
    % lowpass filter
    f = fft(z);
    n = round(4000/continuous_f/2*L);
    f(n+2:L-n) = 0;
    z = real(ifft(f));

end
function y = increaseSampleRate( signal , f1 , f2) % change sample rate from f1 to f2
    L = length(signal);
    T = L/f1 ; % duration of signal
    Ly = round(T*f2) ; % length of y
    t = round((0:1:L)*f2/f1)+1; % scale index
    y = zeros(Ly,1);
    for i=1:1:L
        y(t(i):t(i+1)-1) = signal(i);
    end
end