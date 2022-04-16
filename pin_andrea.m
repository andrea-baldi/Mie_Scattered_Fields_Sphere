% Function pi_n as described in eq. 4.47 in Bohren and Huffmann.
% Inputs are the order of the vector spherical harmonics, 'n', and the
% cosine of the polar angle, 'cos(theta)'
function p=pin_andrea(n,costheta)
if n==0
    p=0;
elseif n==1
    p=1;
else
    q=[0,1];
    for j=2:n
        q(j+1) = costheta*q(j)*(2*j-1)/(j-1)-j*q(j-1)/(j-1);
    end
    p=q(end);
end
