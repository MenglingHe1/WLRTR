function [G1] = WGGG1(G1,Y1,Z1,P1,P2,U1,s0,sf,par,D)
alpha =  par.alpha;
beta  =  par.beta ;
mu    =  par.mu;
U = zeros(size(G1));
M1 = U;
N1 = U;
V1 = U;
iter = 0;
while iter < 10
        iter = iter + 1; 
        % R1
        J1 = double(ttm(tensor(G1),D,2)) - (M1./mu);
 R1 = sign(J1).* max(abs(J1) - alpha.*(1./(abs(J1)+eps))./ mu , 0);
  %  R1 =  soft( J1 , (tau/beta)*(1./(abs(J1)+eps)) );
   
 [ G1] = WCG1(G1,Y1,Z1,P1,P2,U1,R1,M1,par,s0,sf,D,V1,N1); 
 [ V1] =  my_wtnn2( G1 - (N1/mu) , beta/mu  );

      M1 = M1 + mu*(R1 - double(ttm(tensor(G1),D,2)));
        N1 = N1 + mu * (V1 - G1);
        
 
end


