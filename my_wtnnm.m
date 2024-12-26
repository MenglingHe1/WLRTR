function  [X] =  my_wtnnm( Y, tua,C )
    [U,SigmaY,V] = svd(full(Y),'econ');
    SigmaY=diag(SigmaY);  
    X=zeros(size(Y));

    TempC  = tua*C./( SigmaY + eps );% Weight vector
    svp = length(find(SigmaY>TempC));
    if svp>=1
        SigmaX = SigmaY(1:svp)-TempC(1:svp);
    X =  U(:,1:svp)*diag(SigmaX)*V(:,1:svp)';   
    end
return;
