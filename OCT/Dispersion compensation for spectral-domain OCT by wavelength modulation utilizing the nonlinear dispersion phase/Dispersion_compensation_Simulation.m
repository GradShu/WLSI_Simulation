clear all; clc;
%Original spectral interference signal
handle  = OCTFileOpen('Default_0193_Mode2D.oct');
disp( OCTFileGetProperty(handle, 'AcquisitionMode') );
disp( OCTFileGetProperty(handle, 'RefractiveIndex') );
disp( OCTFileGetProperty(handle, 'Comment') );
disp( OCTFileGetProperty(handle, 'Study') );
disp( OCTFileGetProperty(handle, 'ExperimentNumber') );
%-----------------------------------------------------------
%reading spectral raw data
NrRawData = OCTFileGetNrRawData(handle);
[RawData, Spectrum] = OCTFileGetRawData(handle, 0);
load('Wavelength.mat','Wavelength');
figure(1);plot(1:2048,Wavelength);
figure(2);plot(Wavelength,Spectrum);
% figure(3);plot(Wavelength,RawData(:,1));
%波数空间采样及线性矫正
%将波长和光谱强度逆序排列
% spectrum = RawData(:,1);
intensity = zeros(1,length(Spectrum));
wavelength = zeros(1,length(Wavelength));
for i = 1:length(Wavelength)
    wavelength(length(wavelength)+1-i) = Wavelength(i);
    intensity(length(intensity)+1-i) = Spectrum(i);
end
%波长映射到波数
sigma = 1./wavelength;
%波数均匀采样
sigma_inter = linspace(sigma(1),sigma(end),length(sigma)); 
%光谱强度三次样条插值
intensity_inter = spline(sigma,intensity,sigma_inter);
figure(4);plot(sigma_inter,intensity_inter);
%光谱干涉信号S(\sigma) = A(\sigma)*cos(4\pi*(z-z0)*\sigma-4\pi*T*n(\sigma)*\sigma)
spec_signal = intensity_inter;
%设置z = 1000um
z = 500;
z0 = 0;
interference_signal = spec_signal.*cos(4*pi*(z-z0).*sigma_inter);%+2*spec_signal;
figure(5);plot(sigma_inter,interference_signal);
%对光谱进行傅里叶变换
Nfft = 2^20;
interference_signal = fft(interference_signal,Nfft);
% interference_signal = fftshift(interference_signal);
deltaSigma = (sigma_inter(end)-sigma_inter(1))/(length(sigma_inter)-1);
deltaZ = 1/(2*Nfft*deltaSigma);
z = (1:Nfft)*deltaZ;
intensity_fft = abs(interference_signal);
angle_fft = angle(interference_signal);
%坐标变换
% z = z-z(end)/2;
figure(6);plot(z,intensity_fft);
figure(7);plot(z,unwrap(angle_fft));

% %取正半部分，对光谱进行傅里叶逆变换。(有问题？)
intensity_fft = intensity_fft(length(intensity_fft)/2+1:end);
% angle_fft = angle_fft(length(angle_fft)/2+1:end);
% z = z(length(z)/2+1:end);
% figure(8);plot(z,intensity_fft);
% figure(9);plot(z,angle_fft);
% intensity_ifft = intensity_fft.*exp(1i*angle_fft);



%对光谱进行傅里叶逆变换
Nifft = 2^20;
intensity_ifft = ifft(intensity_ifft,Nifft);
intensity_ifft = ifftshift(intensity_ifft);
deltaZ = (z(end)-z(1))/(length(z)-1);
deltaSigma = 1/(2*Nifft*deltaZ);
sigma = (1:Nifft)*deltaSigma;
intensity_ifft = abs(intensity_ifft);
angle_ifft = angle(intensity_ifft);
%坐标变换
sigma = sigma-sigma(end)/2;
figure(10);plot(sigma,intensity_ifft);
figure(11);plot(sigma,angle_ifft);