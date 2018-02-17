function output = MM( A, K4 )
%michaelis menten terms
output = zeros(size(A));

for i = 1:size(K4,1)
    var1 = K4(i,1); %
    var2 = K4(i,2);
    output(var1) = K4(i,3)*A(var2)/(K4(i,4) + A(var2));
end
    
end

