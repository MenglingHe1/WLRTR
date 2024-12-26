function [X] = my_wtnn2(Y,tua)

[n1,n2,n3] = size(Y);  
X = zeros(n1,n2,n3);      
Y = fft(Y,[],2);       

n=1;
NSig = 1;
J = min(n1,n2);
C = n*sqrt(J)*NSig^2;

X(:,1,:) = my_wtnnm( squeeze(Y(:,1,:)),tua, C);
% i=2,...,halfn3
halfn2 = round(n2/2);
for i = 2 : halfn2
   X(:,i,:) = my_wtnnm( squeeze(Y(:,i,:)),tua, C);
    X(:,n2+2-i,:) = conj(X(:,i,:));
end

% if n3 is even
if mod(n2,2) == 0
    i = halfn2+1;
    X(:,i,:) = my_wtnnm( squeeze(Y(:,i,:)),tua, C);
end
X = ifft(X,[],2);