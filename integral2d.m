function [I] = integral2d( z ,  h )
%Integrates 2-dimensional array z with grid spacing of dx using 6 point
%rule, edge points are treated specially with no rule.

I = 0;

mesh = size(z);

%middle points
for i = 2:(mesh(1)-1);
    for j = 2:(mesh(2)-1);
        dI = 4/9*z(i,j) + 1/36* (z(i-1,j-1) + z(i-1,j+1) + z(i+1,j-1) + z(i+1,j+1)) +...
            1/9* (z(i+1,j)+z(i-1,j)) + 1/9 * (z(i,j+1)+z(i,j-1));
        
        I = I + dI*h^2;
    end
end

%left and right edge and corner points (much less accurate)

for i = 1:mesh(1);
    dI = z(i,1);
    I = I +dI*h^2;
    dI = z(i,mesh(2));
    I = I + dI*h^2;
end

%top and bottom
for j = 2:(mesh(2)-1);
    dI = z(1,j);
    I = I +dI*h^2;
    dI = z(mesh(1),j);
    I = I +dI*h^2;
end

end



