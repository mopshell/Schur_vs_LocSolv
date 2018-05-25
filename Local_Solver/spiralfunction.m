function  out = spiralfunction(m,n)

A    = 1;
T    = 1;
B    = m;
L    = 1;
R    = n;
dir  = 0;
out  = zeros(m,n);

while (T<=B && L<=R)

  if dir==0
    for i = T:B
      out(i,L) = A;
      A   = A+1;
    end
      L   = L+1;
      dir = 1;
    
  elseif dir==1

    for i = L:R
      out(B,i) = A;
      A   = A+1;
    end
      B   = B-1;
      dir = 2;
    

  elseif dir==2
    for i = B:-1:T
      out(i,R) = A;
      A   = A+1;
    end
      R   = R - 1;
      dir = 3;
    

  elseif dir==3
      for i = R:-1:L
        out(T,i) = A;
        A   = A+1;
      end
        T   = T + 1;
        dir = 0;
      
  end
  
end
end
