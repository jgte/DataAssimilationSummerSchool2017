function [out,extremity_length,M]=e4_diff(y,dt,mode)

  global order

  switch mode
  case 'diff'
    out=diff(y,2,1)/dt/dt;
    M=[1,-2,1];
    %extremity length is only one
    extremity_length=1;
  case 'taylor'
    [out,M]=finite_diff(y,dt,order,2);
    extremity_length=order;
  otherwise
    error(['unknown mode ''',mode,'''.'])
  end

end

function [df,M]=finite_diff(y,dt,n,nderiv)
  %make sure there's enough points to compute the derivative
  assert(nderiv<2*n+1,['Cannot calculate the ',num2str(nderiv),'-th derivative with only order ',num2str(n),'.'])
  %indexes of the data points
  i=transpose(-n:n);
  %values of the power
  p=0:2*n;
  %compute the inverse of the vandermonde matrix
  M=inv(...
    ( i*ones(    1,2*n+1)).^(ones(2*n+1,1)*p).*...
    (dt*ones(2*n+1,2*n+1)).^(ones(2*n+1,1)*p)...
  ).*(...
    transpose(factorial(0:2*n))*ones(1,2*n+1)...
  );
  %pick the row corresponding to the requested derivative and apply sign change (when relevant)
  M=M(nderiv+1,:)*(-1)^nderiv; 
  %make room for outputs
  df=zeros(size(y,1)-2*n,size(y,2));
  %convolute the filter coefficients over the data, one column at a time
  for i=1:size(y,2)
    df(:,i)=conv(y(:,i),M,'valid');
  end
end

