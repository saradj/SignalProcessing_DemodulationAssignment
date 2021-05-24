close all
clear


%---- Q1 and Q2

filename = 'sig23b_lowpass500hz_v2.mp4';
[y, Fs] = audioread(filename); %reading the signal
N = length(y);

%-- uncomment to play the sound!
%{
player = audioplayer(y,Fs);
play(player);
%}

%---- Q3
Y = fft(y, N); 

%---- Q4
figure
plot(abs(Y));
title('magnitude of fft');

%--------getting the recovered music 1

impulse = zeros(1, N);
f0 = 255000;
impulse(f0) = 1; %impulse at the middle of the first music
impulse(N - f0 + 1) = 1; % a symetric impulse for the negative side
impulse = reshape(impulse, N, 1); % needed to be able to do the convolution!

Shifted = cconv(Y,impulse, N); %shifting the signal by a circular convolution with the impulse train

figure
plot(abs(Shifted));
title('shifted music 1');

frequency=(0 : N - 1)*Fs/N; % to get the accurate measure for the frequency needed for finding fmax

figure
plot(frequency(1 : end/2), abs(Y(1 : end/2)));
title('finding max passband frequency using the fft');

fpass = 2250; % found by inspection on the fft plot (half of the width of the "mountain"/piece 1
low_pass1 = lowpass(ifft(Shifted), fpass, Fs); %filtering all the higher frequencies, keeping ony the piece one!

figure
plot(abs(low_pass1));
title('result music 1');

%----- Alternatively, we could have used multiplication with the inverse fft of the impulse and the original signal
%in the time domain and obtain the same result, uncomment to plot and see equivalence

%{
i = ifft(impulse);
mult = i.*y;
mult_res = lowpass(mult, 2250, Fs);
figure
subplot(2,1,1);plot((abs(low_pass1)));
title('result after convolution in frequency domain');
subplot(2,1,2);plot((abs(mult_res)));
title('result after multiplication in time domain');
print -djpeg compare.jpg
%}

audiowrite('music1_recovered.mp4' , real(low_pass1) , Fs ); %using only the real component since the imaginary one is negligible!

pause(10)
%--------getting the recovered music 2

impulse2 = zeros(1, N);
f1 = 680000;
impulse2(f1) = 1; %impulse at the middle of the first music
impulse2(N - f1 + 1) = 1; % a symetric impulse for the negative side
impulse2 = reshape(impulse2, N, 1); % needed to be able to do the convolution!

Shifted2 = cconv(Y,impulse2, N); %shifting the signal by a circular convolution with the impulse train

figure
plot(abs(Shifted2));
title('shifted music 2');

fpass2 = 2643; % found by inspection on the fft plot (half of the width of the "mountain"/piece 2
low_pass2 = lowpass(ifft(Shifted2), fpass2, Fs); %filtering all the higher frequencies, keeping ony the piece two!

figure
plot(abs(low_pass2));
title('result music 2');

audiowrite('music2_recovered.mp4' , real(low_pass2) , Fs ); %using only the real component since the imaginary one is negligible!

pause(10)
%---- Q5
[recovered1, Fs_rcv1] = audioread('music1_recovered.mp4');
[recovered2, Fs_rcv2] = audioread('music2_recovered.mp4');

%--- getting the magnitude and phase of the fft of recovered music 1 
fft_1 = fft(recovered1);
Magnitude1 = abs(fft_1);
Phase1 = angle(fft_1);

%--- getting the magnitude and phase of the recovered music 1
fft_2 = fft(recovered2);
Magnitude2 = abs(fft_2);
Phase2 = angle(fft_2);

%---synthesizing a signal whose fft has same magnitude as music 1 and same
% phase as music 2
M1_A2 = Magnitude1.*exp(1i*Phase2);
new_music1 = ifft(M1_A2);
audiowrite('newMusic1.mp4', real(new_music1) , Fs_rcv1 );

pause(10)

%---synthesizing a signal whose fft has same magnitude as music 2 and same
% phase as music 1
M2_A1 = Magnitude2.*exp(1i*Phase1);
new_music2 = ifft(M2_A1);
audiowrite('newMusic2.mp4', real(new_music2) , Fs_rcv2 );
