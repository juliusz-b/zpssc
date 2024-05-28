function [xcov_out] = calcXcov(data_in, data_original, should_ref, xcov_block)

if nargin<4
    xcov_block{1} = 'normal';
    xcov_block{2} = 1;
    xcov_block{3} = 1;
    if nargin<3
        should_ref = false;
        if nargin<2
           data_original = data_in; 
        end
    end
end
if size(data_in,1)==size(data_original,1)
   should_ref=false; 
end


xcov_mode = xcov_block{1};
if ~contains(xcov_mode, {'mean', 'normal'})
    error('Wybrany tryb nie moze byc obsluzony');
end

xcov_samples = xcov_block{2}*xcov_block{3};


MODE = 2;
if MODE ~= 1
   warning('Wykorzystywany jest inny tryb MODE niz 1!') 
end

M = size(data_in, 2);if M<size(data_original, 2);M = size(data_original, 2);end
xcov_out = zeros(size(data_original, 1), M);

if (size(data_in)~=size(data_original))
    if ~should_ref
        error('Rozmiary danych porownywanych i oryginalnych powinny byc takie same!')
    end
end

s_ = size(data_original, 1);

if MODE~=2
    if size(data_original,2)<size(data_in,2)
        data_original = [data_original, zeros(size(data_original,1), size(data_in,2) - size(data_original,2))];
    end
end

for i=1:s_
    
    if i == 5
       xd = 1; 
    end
    
    if should_ref
        if MODE==1
            %aa = xcovBlockCalc(data_in(1,:), data_original(i,:), 500);
            xcov_out(i,:) = myXcov(data_in(1,:), data_original(i,:));
            %xcov_out(i,:) = movingMyXcov(data_in(1,:), data_original(i,:));
        elseif MODE==2
            xcov_temp = xcorr(data_in(1,:), data_original(i,:));
            xcov_out(i,:) = xcov_temp(end/2:end);
        else
            xcov_temp = conv(data_in(1,:)-mean(data_in(1,:)), flip(data_original(i,:)-mean(data_original(i,:))));
            xcov_out(i,:) = xcov_temp(end/2:end);
        end
    else
        if MODE==1
            xcov_out(i,:) = myXcov(data_in(i,:), data_original(i,:));
            %xcov_out(i,:) = movingMyXcov(data_in(i,:), data_original(i,:));
        elseif MODE==2
            xcov_temp = xcorr(data_in(i,:), data_original(i,:));%/length(data_original(i,:))*2;
            xcov_out(i,:) = xcov_temp(end/2:end);
        else
            xcov_temp = conv(data_in(1,:)-mean(data_in(1,:)), flip(data_original(i,:)-mean(data_original(i,:))));
            xcov_out(i,:) = xcov_temp(end/2:end);
        end
    end
end

end

