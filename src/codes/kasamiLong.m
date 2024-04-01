function [out] = kasamiLong(m, polyL, initial_polyL)
if m<4||m>23&&(m~=21||m~=22);error('Prbs class unknown');end

if mod(m,4)~=2;error('Wrong power of Long Kasami code!');end

L = 2^m-1;

outL = lfsr(polyL, initial_polyL, L);
q = 2^(m/2)+1;
d_out = repmat(outL,1,q);
d_out = d_out(q:q:end);

qPrim = 2^(m/2+1)+1;
v_out = repmat(outL,1,qPrim);
v_out = v_out(qPrim:qPrim:end);

%u - outL
%w = d_out
%v = v_out


K = 2^m;
M = 2^(m/2)-1;
it = 0;
for k=0:K
    for mM=0:M
        it = it + 1;
        if k>=0&&k<=2^m-2 && m>=0 && m<= 2^(m/2)-2
            seq = xor(outL, circshift(v_out, k));
            seq = xor(seq, circshift(d_out, mM));
        elseif k==2^m-1 && m>=0 && m<=2^(m/2)-2
            seq = xor(outL, circshift(d_out, mM));
        elseif k==2^m && m>=0 && m<= 2^(m/2)-2
            seq = xor(v_out, circshift(d_out, mM));
        elseif k>=0 && k<= 2^m-2 && m==2^(m/2)-1
            seq = xor(outL, circshift(v_out, k));
        elseif k==2^m-1 && m ==2^(m/2)-1
            seq = xor(v_out);
        elseif k==2^m && m==2^(m/2)-1
            seq = xor(d_out);
        else
            error('You should be here!');
        end
        out(L*(it-1)+1:it*L) = double(seq);
    end
end

end

