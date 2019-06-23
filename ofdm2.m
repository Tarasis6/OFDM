clear all
close all
load /home/taras/tpi/ofdm/Demchenko.mat

T_el = 1008;
Fs = 2048/252;
N = Fs*T_el;
GI = 1/4 * N;
Ns = N + GI


a = rfeDump(N: 19*Ns+N-1 +GI  );
b = rfeDump(1 : 19*Ns +GI );
c = a .* conj(b);
n = 0;
while n < length(c) - GI 
    triangle(n+1) = sum( c(1 + n : n + GI )) ;
    n = n + 1; 
end
figure
plot(abs(triangle));

triangle_shape = reshape(triangle, Ns, [] );

figure
plot(abs(triangle_shape()));

triangle_sum = sum(triangle_shape,2);
[y,x] = max(abs(triangle_sum));
figure
plot(abs(triangle_sum));

dfi = angle(triangle_sum(x))
n = 0;
while n < length(rfeDump)
    rfeDump_cor(n+1) = rfeDump(n+1)*exp(-1i*(n)*dfi/N); 
    n = n + 1; 
end
a_cor = rfeDump_cor(N: 19*Ns+N-1 +GI  );
b_cor = rfeDump_cor(1 : 19*Ns +GI );
c_cor = a_cor .* conj(b_cor);
n= 0;
while n < length(c_cor) - GI 
    triangle1(n+1) = sum( c_cor(1 + n : n + GI )) ;
    n = n + 1; 
end

triangle_shape1 = reshape(triangle1, Ns, [] );

triangle_sum1 = sum(triangle_shape1,2);
dfi2 = angle(triangle_sum1(x))

 symbols = reshape(rfeDump_cor(x:x+24*Ns-1),Ns,[]);
 symbols = symbols(GI+1:end,:);
 S_symbols = fftshift(fft(symbols),1);
  %S_symbols = S_symbols(((N-5616)/2 : (N-5616)/2+5616-1),: );
% %  S = fftshift(fft(symbols(5,:)));
% %  ffff = rfeDump_cor(x :  x+N-1);
% %  S = fftshift(fft(rfeDump_cor(x :  x+N-1)));
% % ff = [-length(S)/2 : length(S)/2-1]*Fs/length(S);

  figure
  plot(abs(S_symbols(:,9)))


 pnSequence = comm.PNSequence('Polynomial',[11 2 0], ...
    'InitialConditions',[1 1 1 1 1 1 1 1 1 1 1],'SamplesPerFrame',5616);

h = -step(pnSequence)*2+1;
h = h';

m1 = [1,0,0,0,0,0,0,0,0,0,0,0];
m2 = [0,0,0,1,0,0,0,0,0,0,0,0];
m3 = [0,0,0,0,0,0,1,0,0,0,0,0];
m4 = [0,0,0,0,0,0,0,0,0,1,0,0];

m1 = repmat(m1,1,468);
m2 = repmat(m2,1,468);
m3 = repmat(m3,1,468);
m4 = repmat(m4,1,468);

mask1 = m1.*h;
mask2 = m2.*h;
mask3 = m3.*h;
mask4 = m4.*h;


mask1 = [zeros(1,1288) mask1 zeros(1,1288)];
mask2 = [zeros(1,1288) mask2 zeros(1,1288)];
mask3 = [zeros(1,1288) mask3 zeros(1,1288)];
mask4 = [zeros(1,1288) mask4 zeros(1,1288)];

M1 = (fft(mask1'));
M2 = (fft(mask2'));
M3 = (fft(mask3'));
M4 = (fft(mask4'));

cor1 = sum(abs((ifft(fft(S_symbols(:,1:4:end)).* conj(M1)))),2);
cor2 = sum(abs((ifft(fft(S_symbols(:,1:4:end)).* conj(M2)))),2);
cor3 = sum(abs((ifft(fft(S_symbols(:,1:4:end)).* conj(M3)))),2);
cor4 = sum(abs((ifft(fft(S_symbols(:,1:4:end)).* conj(M4)))),2);


figure
subplot(4,1,1); plot((cor1)); title('Mask1')
subplot(4,1,2); plot((cor2)); title('Mask2')
subplot(4,1,3); plot((cor3)); title('Mask3')
subplot(4,1,4); plot((cor4)); title('Mask4')

[y1,x1] = max(cor1);
[y2,x2] = max(cor2);
[y3,x3] = max(cor3);
[y4,x4] = max(cor4);



S_symbols = S_symbols((N-5616)/2+x1:(N-5616)/2+x1+5616-1, :).';

S_sym = S_symbols(1,:) .* mask1;

Xp = mask1.*PS*4/3;
Yp = S_sym(1,1:12:end);
hp = Yp./Xp(1:12:end);
x_inter = 1:12:5616;
xq = 1:5616;    
iterp_pilot = interp1(x_inter, hp, xq, 'spline');
subplot(1,2,1)
plot( abs(iterp_pilot))
SSS  = S_symbols(1,:)./iterp_pilot;
subplot(1,2,2)
plot( real(SSS),imag(SSS),'.')
axis([-3,3,-3,3]);
        










        