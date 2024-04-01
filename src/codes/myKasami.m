function [out] = myKasami(m, type)
if isempty(find([4,6,8,10,12,14,16,18]==m, 1, 'first'));error('Wrong Kasami m number!');end
if nargin<2
   type = 'short'; 
end

if ~contains({'long', 'short'}, type);error('Wrong Kasami type!');end

switch m
    case 4
        poly1 = oct2poly(23);
    case 6
        poly1 = oct2poly(103);
    case 8
        poly1 = oct2poly(435);
    case 10
        poly1 = oct2poly(2011);
    case 12
        poly1 = oct2poly(10123);
    case 14
        poly1 = oct2poly(42103);
    case 16
        poly1 = oct2poly(210013);
    case 18
        poly1 = oct2poly(1000201);
end

switch type
    case 'long'
    temp = kasamiLong(m, poly1, [zeros(1,length(poly1)-2) 1]);
    M = (2^(m/2))*(2^(m)+1);
    case 'short'
    temp = kasamiShort(m, poly1, [zeros(1,length(poly1)-2) 1]);
    M = 2^(m/2);
end


n = 2^m-1;
out = [];
for i=1:M
    out(i, :) = temp(n*(i-1)+1:n*i);
end

% figure();
% for i=1:M
%     xcv = crossCov(out(3,:), out(i,:));
%     plot(xcv, 'o');
%     unique(xcv)
%     
%     pause(0.05)
% end

end

