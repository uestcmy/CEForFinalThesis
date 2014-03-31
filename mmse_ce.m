function[H_MMSE]=mmse_ce(H_est,Nfft,Nps,h,SNR)
    k = 0:length(h)-1;
    hh = h*h';
    tmp = h(1,:).*conj(h(1,:)).*k;
    r = sum(tmp)/hh;
    r2 = tmp*k.'/hh;
    tao_rms = sqrt(r2-r*r);
    df = 1/Nfft;
    j2pi_tao_df = i*2*pi*tao_rms*df;
    K1 = repmat([0:Nfft-1].',1,Nfft);
    K2 = repmat([0:Nfft-1],Nfft,1);
    rf = 1./(1+j2pi_tao_df*(K1-K2));
    K3 = repmat([0:Nfft-1].',1,Nfft);
    K4 = repmat([0:Nfft-1],Nfft,1);
    rf2 = 1./(1+j2pi_tao_df*(K3-K4).*Nps);
    Rpp = rf2+eye(length(H_est),length(H_est))/SNR;
    H_MMSE = rf/(Rpp)*H_est;
    H = fft(h,Nfft);
    figure();subplot(2,1,1);plot(10*log10(abs(H.*conj(H)))),hold on,stem(10*log10(H_est.*conj(H_est))),plot(10*log10(H_MMSE.*conj(H_MMSE)),'r:+');
    title(['CE for SISO tao_rms']);legend('True Channel','LS','MMSE',4);
    xlabel( ['tao_rms = ' num2str(tao_rms)]);
    subplot(2,1,2);stem(abs(h));axis([-0.2 4 -0.5 1]);
end