
clc
clear  

addpath( 'newdata','quality_assess', 'TR_function', 'Toolbox')
load('lemon slices.mat')
S=HSI; 
 %F=F;
S=double(S);
S=S/max(S(:));
sf  =  8 ;
s0= sf/2 ;

  [M, N, L]=size(S); 
  S1=hyperConvert2D(S);
F = [2  1  1  1  1  1  0  0  0  0  0  0  0  0  0  0  2  6 11 17 21 22 21 20 20 19 19 18 18 17 17;...
     1  1  1  1  1  1  2  4  6  8 11 16 19 21 20 18 16 14 11  7  5  3  2  2  1  1  2  2  2  2  2;...
     7 10 15 19 25 29 30 29 27 22 16  9  2  0  0  0  0  0  0  0  1  1  1  1  1  1  1  1  1  1  1];
  for band = 1:3
        div = sum(F(band,:));
        for i = 1:31
            F(band,i) = F(band,i)/div;
        end
  end

 psf=fspecial('gaussian',[7,1],2);
 % psf=ones(8,1)/8;
  
 BW1=psf2otf(psf,[M 1]);
 S_w=ifft(fft(S).*repmat(BW1,1,N,L)); 
 BH1=psf2otf(psf,[N 1]);
aa=fft(permute(S_w,[2 1 3]));
  S_h=(aa.*repmat(BH1,1,M,L));
 S_h= permute(ifft(S_h),[2 1 3]);  
  Y_h=S_h(s0:sf:end,s0:sf:end,:);
  Y_h_bar=hyperConvert2D(Y_h);
  
  SNRh=  25 ;
sigmam = sqrt(sum(Y_h_bar(:).^2)/(10^(SNRh/10))/numel(Y_h_bar));
rng(10,'twister')
   Y_h_bar = Y_h_bar+ sigmam*randn(size(Y_h_bar));
HSI=hyperConvert3D(Y_h_bar,M/sf, N/sf );
%save('pavia_HSI.mat')

 %%  simulate HR-MSI
 rng(10,'twister')
Y = F*S1;
SNRm=  30  ;
sigmam = sqrt(sum(Y(:).^2)/(10^(SNRm/10))/numel(Y));
Y = Y+ sigmam*randn(size(Y));
 MSI=hyperConvert3D(Y,M,N);
 %save('pavia_MSI.mat')

%% FSTRD_TV_PAM
%%%%%%%%%%%%%%%%%%%%%%   参数设置     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 45
r    =   5 ;
par.TRrank = [5,300,5];
par.lambda =  1.3; 
par.alpha  =   0.0001;
par.beta   =    0.8 ;
par.mu     =   0.03; 
par.rho    =  1.1 ;

%%%%%%%%%%%%%%%%%%%%%%   模型算法   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Z0,runtime] = WLRTR( HSI,MSI,F,BW1,BH1,S,par ,s0,sf);
[psnr4,rmse4,ergas4,sam4,uiqi4,ssim4,DD4,CC4] = quality_assessment(Z0*255, S*255, 0, 1.0/sf);
fprintf('================== Result =====================\n');
fprintf(' %8.8s  %5.4s  %5.4s  %6.5s  %4.3s  %5.4s  \n','method','PSNR', 'SSIM', 'ERGAS', 'SAM', 'UIQI' );
fprintf(' %8.8s  %2.3f  %2.3f  %2.3f  %2.3f  %2.3f  \n','Ours', psnr4, ssim4,ergas4,sam4,uiqi4);
