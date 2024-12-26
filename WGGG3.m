function [G3] = WGGG3(G3,Y3,Z3,T1,T2,U3,par,D)
alpha =  par.alpha;
beta  =  par.beta ;
mu    =  par.mu;
U = zeros(size(G3));
M3 = U;
N3 = U;
V3 = U;
iter = 0;
while iter < 10
        iter = iter + 1; 
        % R1
        J3 = double(ttm(tensor(G3),D,2)) - (M3/mu);
   R3 = sign(J3).* max(abs(J3) - alpha.*(1./(abs(J3)+eps))./ mu , 0);
%  R3 =  soft( J3 , (tau/beta)*(1./(abs(J3)+eps)) );
  [G3] = WCG3(G3,Y3,Z3,T1,T2,U3,R3,M3,par,D,V3,N3);
 [ V3] =  my_wtnn2(   G3- (N3/mu) , beta/mu  );
  M3 = M3 + mu*(R3 - double(ttm(tensor(G3),D,2)));
        N3 = N3 + mu * (V3 - G3);
    
end

