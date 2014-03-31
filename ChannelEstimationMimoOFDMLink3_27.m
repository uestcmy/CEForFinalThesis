clc,clear all,close all
%% Introduction
% MIMO-OFDM仿真链路，发送端单天线，收端单用户四天线，均采用32个子载波的OFDM传输
% 用户信道均为4个taps的时不变信道，经过归一处理
% 干扰信道均为5个taps的时不变信道，经过归一处理

%% Initialization
% Antenna init**********
Nrx = 4;
Ntx = 1;
% OFDM Init**********
Nsym = 3;
Nfft = 32;
Ng = 6;
Nofdm = Nfft + Ng;
% Modulator Init**********
M = 16;
Mod16QAMObject = comm.RectangularQAMModulator(M);
Mod16QAMObject.NormalizationMethod = 'Average power';
Mod16QAMObject.AveragePower = 1;
ModQPSKObject = comm.QPSKModulator;
% Channel Init**********
SIR = 15;
SNR = 15;
UserChannelTaps = 3;
UserChannel = randn(Nrx,UserChannelTaps) + 1i*randn(Nrx,UserChannelTaps);
for k1 = 1:1:Nrx
    UserChannel(k1,:) = UserChannel(k1,:)/sqrt(abs(UserChannel(k1,:)*UserChannel(k1,:)'));
end
InterChannelTaps = 2;
InterChannel = randn(Nrx,InterChannelTaps) + 1i*randn(Nrx,InterChannelTaps);
for k1 = 1:1:Nrx
    InterChannel(k1,:) = InterChannel(k1,:)/sqrt(abs(InterChannel(k1,:)*InterChannel(k1,:)'));
end

%% Transmitter
% User BS**********  baseband signal
OriData = randi([0,M-1],Nsym*Nfft,1);  %uniform distribution[0,15]  100*32
ModData = step(Mod16QAMObject,OriData);
ModDataRE = reshape(ModData,Nfft,Nsym);
TimeDomDataRE = ifft(ModDataRE,Nfft)*sqrt(Nfft);
TimeDomDataRE_WithCp = [TimeDomDataRE(end-Ng+1:end,:);TimeDomDataRE];
UserTransData = reshape(TimeDomDataRE_WithCp,Nofdm*Nsym,1);
% Interference BS**********
InterOriData = randi([0,3],Nsym*Nfft,1);
InterModData = step(ModQPSKObject,InterOriData);
InterModDataRE = reshape(InterModData,Nfft,Nsym);
InterTimeDomDataRE = ifft(InterModDataRE,Nfft)*sqrt(Nfft);
InterTimeDomDataRE_WithCp = [InterTimeDomDataRE(end-Ng+1:end,:);InterTimeDomDataRE];
InterTransData = reshape(InterTimeDomDataRE_WithCp,Nofdm*Nsym,1);
% Power Allocation
InterPower = 10^(-SIR/10);
NoisePower = 10^(-SNR/10);

%% Channel
RecDataDomin = zeros(Nrx,Nofdm*Nsym);4
for k1 = 1:1:Nrx
    RecDataUserPart = filter(UserChannel(k1,:),1,UserTransData);
    RecDataInterPart = sqrt(InterPower)*filter(InterChannel(k1,:),1,InterTransData);
    RecDataNoisePart = sqrt(NoisePower/2)*(randn(length(UserTransData),1)+1i*randn(length(UserTransData),1));
    RecDataDomin(k1,:) = (RecDataUserPart + RecDataInterPart + RecDataNoisePart).';
end

%% Receiver
% ReMove OFDM Operation
for k1 = 1:1:Nrx
    OperationTimeDomData = RecDataDomin(k1,:);
    OperationTimeDomDataRE_WithCp = reshape(OperationTimeDomData,Nofdm,Nsym);
    OperationDataTimeDomRE = OperationTimeDomDataRE_WithCp(Ng+1:Nofdm,:);
    OperationDataRE(k1,:,:) = fft(OperationDataTimeDomRE,Nfft)/sqrt(Nfft);
end
% ChannelEstimation method 1:LS

for k1 = 1:1:Nrx
    YpilotSymbol = squeeze(OperationDataRE(k1,:,:));
    H_est = mean(YpilotSymbol./ModDataRE,2);
    H = fft(UserChannel(k1,:),Nfft);
    H_mmse = mmse_ce(H_est,Nfft,1,UserChannel(k1,:),SNR);

end
% % ChannelEstimation method 2:MMSE
% k = 0:UserChannelTaps-1;
% hh = UserChannel(1,:)*UserChannel(1,:)';
% tmp = UserChannel(1,:).*conj(UserChannel(1,:)).*k;
% r = sum(tmp)/hh;
% r2 = tmp*k.'/hh;
% tao_rms = sqrt(r2-r*r);
% df = 1/Nfft;
% j2pi_tao_df = i*2*pi*tao_rms*df;
% K1 = repmat([0:Nfft-1].',1,Nfft);
% K2 = repmat([0:Nfft-1],Nfft,1);
% rf = 1./(1+j2pi_tao_df*(K1-K2));
% K3 = repmat([0:Nfft-1].',1,Nfft);
% K4 = repmat([0:Nfft-1],Nfft,1);
% rf2 = 1./(1+j2pi_tao_df*(K3-K4));
% Rpp = rf2+eye(length(H_est),length(H_est))/SNR
% 
% H_mmse = rf*inv(Rpp)*H_est;
% figure(), plot(10*log10(abs(H.*conj(H)))),hold on,title('test mmse');plot(10*log10(H_mmse.*conj(H_mmse)),'r:+');plot(10*log10(H_est.*conj(H_est)),'b:o');
%stem([1:32],[10*log10(H_est.*conj(H_est));10*log10(H_est.*conj(H_est))]);title('test mmse');
%figure(),plot(10*log10(abs(H.*conj(H)))),hold on,st = stem(10*log10(H_mmse.*conj(H_mmse))),set(st,'MarkerFaceColor','blue')
% ChannelEstimation method 3:Signal Estimation Technic

% Symbol DeModulator And Dectection

% BER Or MSE Cal