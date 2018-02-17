function output = NL( A, K3 )
%nonlinear (quadratic) terms
output = zeros(size(A));

for i = 1:size(K3,1)
    var1 = K3(i,1); %
    var2 = K3(i,2);
    var3 = K3(i,3);
    output(var1) = K3(i,4)*A(var2)*A(var3);
end
    
end

