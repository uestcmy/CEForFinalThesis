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
Nsym = 10000;
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
% User BS**********
OriData = randi([0,M-1],Nsym*Nfft,1);
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
RecDataDomin = zeros(Nrx,Nofdm*Nsym);
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
    figure(k1),plot(10*log10(abs(H.*conj(H)))),hold on,stem(10*log10(H_est.*conj(H_est)))
end
% ChannelEstimation method 2:MMSE 

% ChannelEstimation method 3:Signal Estimation Technic

% Symbol DeModulator And Dectection

% BER Or MSE Cal






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test 1:The Likelihood Function
% H_Test = fft(UserChannel,Nfft,2);
% G_Test = fft(InterChannel,Nfft,2);
% %Part I:
% LikelihoodPart1 = -size(OperationDataRE,1)*size(OperationDataRE,2)*size(OperationDataRE,3)*log(sqrt(NoisePower));
% %Part II:
% LikelihoodPart2 = 0;
% for k1 = 1:1:size(OperationDataRE,2)
%     LikelihoodPart2 = LikelihoodPart2 - size(OperationDataRE,3)*log(1+G_Test(:,k1)'*G_Test(:,k1)/NoisePower);
% end
% %Part III:
% LikelihoodPart3 = 0;
% for k1 = 1:1:size(OperationDataRE,2)
%     for k2 = 1:1:size(OperationDataRE,3)
%         ErrEveryRE = OperationDataRE(:,k1,k2) - H_Test(:,k1)*ModDataRE(k1,k2);
%         LikelihoodPart3 = LikelihoodPart3 -.....
%             (1/NoisePower)*ErrEveryRE'*inv(eye(Nrx)+G_Test(:,k1)*G_Test(:,k1)'/NoisePower)*ErrEveryRE;
%     end
% end
% Likelihood = real(LikelihoodPart1 + LikelihoodPart2 + LikelihoodPart3);
% %Random Experience
% for ExperTime = 1:1:1000
%     fprintf('Current Experience %d.\n',ExperTime);
%     RandUserChannelTaps = 3;
%     RandUserChannel = randn(Nrx,RandUserChannelTaps) + 1i*randn(Nrx,RandUserChannelTaps);
%     for k1 = 1:1:Nrx
%         RandUserChannel(k1,:) = RandUserChannel(k1,:)/sqrt(abs(RandUserChannel(k1,:)*RandUserChannel(k1,:)'));
%     end
%     RandInterChannelTaps = 2;
%     RandInterChannel = randn(Nrx,RandInterChannelTaps) + 1i*randn(Nrx,RandInterChannelTaps);
%     for k1 = 1:1:Nrx
%         RandInterChannel(k1,:) = RandInterChannel(k1,:)/sqrt(abs(RandInterChannel(k1,:)*RandInterChannel(k1,:)'));
%     end
%     RandNoisePower = randn^2;
%     RandH_Test = fft(RandUserChannel,Nfft,2);
%     RandG_Test = fft(RandInterChannel,Nfft,2);
%     %Part I:
%     RandLikelihoodPart1 = -size(OperationDataRE,1)*size(OperationDataRE,2)*size(OperationDataRE,3)*log(sqrt(RandNoisePower));
%     %Part II:
%     RandLikelihoodPart2 = 0;
%     for k1 = 1:1:size(OperationDataRE,2)
%         RandLikelihoodPart2 = RandLikelihoodPart2 - .....
%             size(OperationDataRE,3)*log(1+RandG_Test(:,k1)'*RandG_Test(:,k1)/RandNoisePower);
%     end
%     %Part III:
%     RandLikelihoodPart3 = 0;
%     for k1 = 1:1:size(OperationDataRE,2)
%         for k2 = 1:1:size(OperationDataRE,3)
%             RandErrEveryRE = OperationDataRE(:,k1,k2) - RandH_Test(:,k1)*ModDataRE(k1,k2);
%             RandLikelihoodPart3 = RandLikelihoodPart3 -.....
%                 (1/RandNoisePower)*RandErrEveryRE'*inv(eye(Nrx)+RandG_Test(:,k1)*RandG_Test(:,k1)'/RandNoisePower)*RandErrEveryRE;
%         end
%     end
%     RandLikelihood(ExperTime) = real(RandLikelihoodPart1 + RandLikelihoodPart2 + RandLikelihoodPart3);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test2 : Improved Max Likelihood Function
% Eita = reshape(InterChannel,Nrx*InterChannelTaps,1)/sqrt(NoisePower);
% Seta = reshape(UserChannel,Nrx*UserChannelTaps,1);
% FFTMatrix = zeros(Nfft,Nfft);
% for k1 = 1:1:Nfft
%     for k2 = 1:1:Nfft
%         FFTMatrix(k1,k2) = exp(-1i*2*pi*(k1-1)*(k2-1)/Nfft);
%     end
% end
% %ImproLikelihoodPart I:
% ImproLikelihoodPart1 = -size(OperationDataRE,1)*size(OperationDataRE,2)*size(OperationDataRE,3)*log(sqrt(NoisePower));
% %ImproLikelihoodPart II:
% ImproLikelihoodPart2 = 0;
% for k1 = 1:1:size(OperationDataRE,2)
%     UGaink1 = kron(FFTMatrix(k1,1:InterChannelTaps),eye(Nrx,Nrx));
%     ImproLikelihoodPart2 = ImproLikelihoodPart2 - .....
%         size(OperationDataRE,3)*log(1+(UGaink1*Eita)'*(UGaink1*Eita));
% end
% %ImproLikelihoodPartIII:
% ImproLikelihoodPart3 = 0;
% for k1 = 1:1:size(OperationDataRE,2)
%     for k2 = 1:1:size(OperationDataRE,3)
%         UGaink1 = kron(FFTMatrix(k1,1:InterChannelTaps),eye(Nrx,Nrx));
%         HGaink1 = kron(FFTMatrix(k1,1:UserChannelTaps),eye(Nrx,Nrx));
%         ErrEveryRE = OperationDataRE(:,k1,k2) - ModDataRE(k1,k2)*HGaink1*Seta;
%         ImproLikelihoodPart3 = ImproLikelihoodPart3 +.....
%             (1/NoisePower)*ErrEveryRE'*(UGaink1*Eita)*(UGaink1*Eita)'*ErrEveryRE/.....
%             (1+(UGaink1*Eita)'*(UGaink1*Eita));
%     end
% end
% %ImproLikelihoodPartIV:
% ImproLikelihoodPart4 = 0;
% for k1 = 1:1:size(OperationDataRE,2)
%     for k2 = 1:1:size(OperationDataRE,3)
%         HGaink1 = kron(FFTMatrix(k1,1:UserChannelTaps),eye(Nrx,Nrx));
%         ErrEveryRE = OperationDataRE(:,k1,k2) - ModDataRE(k1,k2)*HGaink1*Seta;
%         ImproLikelihoodPart4 = ImproLikelihoodPart4 - .....
%             (1/NoisePower)*ErrEveryRE'*ErrEveryRE;
%     end
% end
% ImproLikelihood = ImproLikelihoodPart1 + ImproLikelihoodPart2 + ImproLikelihoodPart3 + ImproLikelihoodPart4;
% %Rand Experience
% for ExperTime = 1:1:2000
%     fprintf('Current Experience %d.\n',ExperTime);
%     RandUserChannelTaps = 3;
%     RandUserChannel = randn(Nrx,RandUserChannelTaps) + 1i*randn(Nrx,RandUserChannelTaps);
%     for k1 = 1:1:Nrx
%         RandUserChannel(k1,:) = RandUserChannel(k1,:)/sqrt(abs(RandUserChannel(k1,:)*RandUserChannel(k1,:)'));
%     end
%     RandInterChannelTaps = 2;
%     RandInterChannel = randn(Nrx,RandInterChannelTaps) + 1i*randn(Nrx,RandInterChannelTaps);
%     for k1 = 1:1:Nrx
%         RandInterChannel(k1,:) = RandInterChannel(k1,:)/sqrt(abs(RandInterChannel(k1,:)*RandInterChannel(k1,:)'));
%     end
%     RandNoisePower = randn^2;
%     RandEita = reshape(RandInterChannel,Nrx*RandInterChannelTaps,1)/sqrt(RandNoisePower);
%     RandSeta = reshape(RandUserChannel,Nrx*RandUserChannelTaps,1);
%     FFTMatrix = zeros(Nfft,Nfft);
%     for k1 = 1:1:Nfft
%         for k2 = 1:1:Nfft
%             FFTMatrix(k1,k2) = exp(-1i*2*pi*(k1-1)*(k2-1)/Nfft);
%         end
%     end
%     %RandImproLikelihoodPart I:
%     RandImproLikelihoodPart1 = -size(OperationDataRE,1)*size(OperationDataRE,2)*size(OperationDataRE,3)*log(sqrt(RandNoisePower));
%     %RandImproLikelihoodPart II:    
%     RandImproLikelihoodPart2 = 0;
%     for k1 = 1:1:size(OperationDataRE,2)
%         RandUGaink1 = kron(FFTMatrix(k1,1:RandInterChannelTaps),eye(Nrx,Nrx));
%         RandImproLikelihoodPart2 = RandImproLikelihoodPart2 - .....
%         size(OperationDataRE,3)*log(1+(RandUGaink1*RandEita)'*(RandUGaink1*RandEita));
%     end
%     %RandImproLikelihoodPartIII:
%     RandImproLikelihoodPart3 = 0;
%     for k1 = 1:1:size(OperationDataRE,2)
%         for k2 = 1:1:size(OperationDataRE,3)
%             RandUGaink1 = kron(FFTMatrix(k1,1:RandInterChannelTaps),eye(Nrx,Nrx));
%             RandHGaink1 = kron(FFTMatrix(k1,1:RandUserChannelTaps),eye(Nrx,Nrx));
%             RandErrEveryRE = OperationDataRE(:,k1,k2) - ModDataRE(k1,k2)*RandHGaink1*RandSeta;
%             RandImproLikelihoodPart3 = RandImproLikelihoodPart3 +.....
%             (1/RandNoisePower)*RandErrEveryRE'*(RandUGaink1*RandEita)*(RandUGaink1*RandEita)'*RandErrEveryRE/.....
%             (1+(RandUGaink1*RandEita)'*(RandUGaink1*RandEita));
%         end
%     end
%     %RandImproLikelihoodPartIV:
%     RandImproLikelihoodPart4 = 0;
%     for k1 = 1:1:size(OperationDataRE,2)
%         for k2 = 1:1:size(OperationDataRE,3)
%             RandHGaink1 = kron(FFTMatrix(k1,1:RandUserChannelTaps),eye(Nrx,Nrx));
%             RandErrEveryRE = OperationDataRE(:,k1,k2) - ModDataRE(k1,k2)*RandHGaink1*RandSeta;
%             RandImproLikelihoodPart4 = RandImproLikelihoodPart4 - .....
%             (1/RandNoisePower)*RandErrEveryRE'*RandErrEveryRE;
%         end
%     end
%     RandImproLikelihood(ExperTime) = real(RandImproLikelihoodPart1 + RandImproLikelihoodPart2 + .....
%         RandImproLikelihoodPart3 + RandImproLikelihoodPart4);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Test 3 : (4.5)
%Test Seta(ML+GLS)
Eita = reshape(InterChannel,Nrx*InterChannelTaps,1);
FFTMatrix = zeros(Nfft,Nfft);
for k1 = 1:1:Nfft
    for k2 = 1:1:Nfft
        FFTMatrix(k1,k2) = exp(-1i*2*pi*(k1-1)*(k2-1)/Nfft);
    end
end
QMatrix = zeros(UserChannelTaps*Nrx,UserChannelTaps*Nrx);
CalQMatrix = zeros(Nrx*UserChannelTaps,1);
for k1 = 1:1:size(OperationDataRE,2)
    UGaink1 = kron(FFTMatrix(k1,1:InterChannelTaps),eye(Nrx,Nrx));
    HGaink1 = kron(FFTMatrix(k1,1:UserChannelTaps),eye(Nrx,Nrx));
    DMinus1k1 = eye(Nrx,Nrx) - (UGaink1*Eita)*(UGaink1*Eita)'/(1+(UGaink1*Eita)'*(UGaink1*Eita));
    sjx = sum(conj(ModDataRE(k1,:)).*ModDataRE(k1,:));
    QMatrix = QMatrix + sjx*HGaink1'*DMinus1k1*HGaink1;
end
for k1 = 1:1:size(OperationDataRE,2)
    UGaink1 = kron(FFTMatrix(k1,1:InterChannelTaps),eye(Nrx,Nrx));
    HGaink1 = kron(FFTMatrix(k1,1:UserChannelTaps),eye(Nrx,Nrx));
    DMinus1k1 = eye(Nrx,Nrx) - (UGaink1*Eita)*(UGaink1*Eita)'/(1+(UGaink1*Eita)'*(UGaink1*Eita));
    sjxy = 0;
    for k2 = 1:1:size(OperationDataRE,3)
        sjxy = sjxy + conj(ModDataRE(k1,k2))*OperationDataRE(:,k1,k2);
    end
    CalQMatrix = CalQMatrix + HGaink1'*DMinus1k1*sjxy;
end
Seta = inv(QMatrix)*CalQMatrix;
%Test NoisePower(ML+GLS)
NoisePowerML = 0;
for k1 = 1:1:size(OperationDataRE,2)
    UGaink1 = kron(FFTMatrix(k1,1:InterChannelTaps),eye(Nrx,Nrx));
    HGaink1 = kron(FFTMatrix(k1,1:UserChannelTaps),eye(Nrx,Nrx));
    DMinus1k1 = eye(Nrx,Nrx) - (UGaink1*Eita)*(UGaink1*Eita)'/(1+(UGaink1*Eita)'*(UGaink1*Eita));
    ErrEveryCarrier = 0;
    for k2 = 1:1:size(OperationDataRE,3)
        ErrEveryCarrier = ErrEveryCarrier + OperationDataRE(:,k1,k2) - ModDataRE(k1,k2)*HGaink1*Seta;
    end
    NoisePowerML = NoisePowerML + (1/size(OperationDataRE,1)/size(OperationDataRE,2)/size(OperationDataRE,3))*.....
        ErrEveryCarrier'*DMinus1k1*ErrEveryCarrier;
end
%Test Eita (ConcenTrated Likelihood Function (4.6))
% Eita = reshape(InterChannel,Nrx*InterChannelTaps,1);
% CalNullSpaceMatrix = zeros(Nrx*InterChannelTaps,1);
Eita = zeros(Nrx*InterChannelTaps,1);
CalNullSpaceMatrix = zeros(Nrx*InterChannelTaps,Nrx*InterChannelTaps);
for k1 = 1:1:size(OperationDataRE,2)
    UGaink1 = kron(FFTMatrix(k1,1:InterChannelTaps),eye(Nrx,Nrx));
    HGaink1 = kron(FFTMatrix(k1,1:UserChannelTaps),eye(Nrx,Nrx));
    Ej = 0;
    for k2 = 1:1:size(OperationDataRE,3)
        eHatjk = OperationDataRE(:,k1,k2) - ModDataRE(k1,k2)*HGaink1*Seta;
        Ej = Ej + 1/(size(OperationDataRE,3))*eHatjk*eHatjk';
    end
    GjEita = ((NoisePowerML*(1+(UGaink1*Eita)'*(UGaink1*Eita))-(UGaink1*Eita)'*Ej*(UGaink1*Eita))*eye(Nrx,Nrx) +.....
        (1+(UGaink1*Eita)'*(UGaink1*Eita))*Ej)/(NoisePowerML*(1+(UGaink1*Eita)'*(UGaink1*Eita))^2);
    CalNullSpaceMatrix = CalNullSpaceMatrix + UGaink1'*GjEita*UGaink1;
end
    


