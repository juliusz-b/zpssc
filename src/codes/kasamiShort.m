function [out] = kasamiShort(m, polyL, initial_polyL)
if m<4||m>23&&(m~=21||m~=22);error('Prbs class unknown');end

if mod(m,2)~=0;error('Wrong power of Kasami code!');end

L = 2^m-1;

outL = lfsr(polyL, initial_polyL, L);
q = 2^(m/2)+1;
d_out = repmat(outL,1,q);
d_out = d_out(q:q:end);

for i=1:2^(m/2)-1
    out(L*(i-1)+1:i*L) = double(xor(outL, circshift(d_out, -(i-1))));
end

out = [outL, out];

end


