function y = array_response_ULA(N,fk,fc,theta)
for m= 0:N-1
    y(m+1) = exp( 1i* pi* fk/fc * ( m*sin(theta) ) );
end
y = y.'/sqrt(N);
end