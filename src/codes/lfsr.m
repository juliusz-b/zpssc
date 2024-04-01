function [out_code] = lfsr(polynomial, beginning_polynomial, number_of_iterations)

if ~isvector(polynomial);error('polynomial must be a vector');end
L = length(polynomial)-1;
if length(beginning_polynomial)~=L;error('beginning_polynomial must have same lengt as polynomial-1');end

reg = beginning_polynomial;

flipped_polynomial = flip(polynomial(1:end-1));
flipped_polynomial_temp = flipped_polynomial;
ones_L = sum(flipped_polynomial);

out_code = zeros(1,number_of_iterations);
for it=1:number_of_iterations
    %for i=1:L
        flipped_polynomial_temp = flipped_polynomial;
        
        if ones_L>1
            for j=1:ones_L-1
                if j==1
                    id1 = find(flip(flipped_polynomial_temp)==1,1,'first');
                    flipped_polynomial_temp = flip(flipped_polynomial_temp);flipped_polynomial_temp(id1)=0;flipped_polynomial_temp = flip(flipped_polynomial_temp);
                    id1 = L-id1+1;
                    id2 = find(flip(flipped_polynomial_temp)==1,1,'first');
                    flipped_polynomial_temp = flip(flipped_polynomial_temp);flipped_polynomial_temp(id2)=0;flipped_polynomial_temp = flip(flipped_polynomial_temp);
                    id2 = L-id2+1;
                    bit = xor(reg(id1), reg(id2));
                else
                    id1 = find(flip(flipped_polynomial_temp)==1,1,'first');
                    flipped_polynomial_temp = flip(flipped_polynomial_temp);flipped_polynomial_temp(id1)=0;flipped_polynomial_temp = flip(flipped_polynomial_temp);
                    id1 = L-id1+1;
                    bit = xor(bit, reg(id1));
                end
            end
            out_code(it) = reg(end);
            reg = circshift(reg, 1);
            reg(1) = bit;
        else
            out_code(it) = reg(end);
            reg = circshift(reg, 1);
            if sum([flipped_polynomial]==1)==0
                reg(1) = reg(end);
            else
                reg(1) = reg([flipped_polynomial]==1);
            end
        end
    %end
end

end

