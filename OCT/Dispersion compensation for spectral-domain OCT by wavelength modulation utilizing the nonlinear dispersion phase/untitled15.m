 sigma_inter1 = sigma+C*(angle_fit-ang_fit)/(-4*pi);
%插值求得处理后的干涉信号,需要修正，首尾坐标选取sigma_inter、sigma_inter1中合适的值
sigma_find = (sigma_inter1>=sigma_inter(1) & sigma_inter1 <=sigma_inter(end));
sigma_inter1 = sigma_inter1(sigma_find);
signal_C = interp1(sigma_inter,interference_signal1,sigma_inter1,'spline');
plot(sigma_inter1,signal_C)
xlabel('Wavenumber');ylabel('Spectral Interference Signal');title("消除色差后");

%对光谱进行傅里叶变换
Nfft = 2^20;
fft_signal = fft(signal_C,Nfft);
interference_signal = fftshift(fft_signal);
deltaSigma = (sigma_inter1(end)-sigma_inter1(1))/(length(sigma_inter1)-1);
deltaZ = 1/(2*Nfft*deltaSigma);
z = (1:Nfft)*deltaZ;
intensity_fft = abs(interference_signal);
angle_fft = angle(interference_signal);

%坐标变换
z = z-z(end)/2;
z = z(length(z)/2+1:end);
intensity_fft = intensity_fft(length(intensity_fft)/2+1:end);
% plot(z,intensity_fft);
threshold = max(intensity_fft)/6;
findpeaks(intensity_fft,z,'MinPeakProminence',threshold,'MinPeakDistance',50,...
 'Annotate','extents');
[pks,locs,w,p] = findpeaks(intensity_fft,z,'MinPeakProminence',threshold,...
                'MinPeakDistance',50,'Annotate','extents');
text(locs,pks+pks/15,num2str((1:numel(pks))'));
xlabel('z');ylabel('Intensity');title("消除色差后");