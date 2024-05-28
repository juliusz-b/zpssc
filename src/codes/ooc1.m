
% Generate a (7,3,1) OOC
codes = generate_ooc(32, 16/4, 6, 1, 8);
figure;plot(xcorr(codes(2,:),codes(3,:)))
hold on;plot(xcorr(codes(2,:)))
%%
% Display the codes
disp('Generated OOC codes:');
disp(codes);

% Visualize the codes
visualize_ooc(codes, 1, 1);

% Generate a (13,3,1) OOC
%codes2 = generate_ooc(13, 3, 1, 1, 2);  % Limit to 2 codes

% Display the codes
%disp('Generated (13,3,1) OOC codes:');
%disp(codes2);

% Visualize the codes
%visualize_ooc(codes2, 1, 1);


function binary_codes = generate_ooc(n, w, lambda_a, lambda_c, max_codes)
    % GENERATE_OOC Generate Optical Orthogonal Codes using greedy algorithm
    % 
    % INPUTS:
    %   n - length of codewords
    %   w - weight of codewords (number of 1's)
    %   lambda_a - autocorrelation constraint
    %   lambda_c - cross-correlation constraint
    %   max_codes - maximum number of codes to generate (optional)
    %
    % OUTPUT:
    %   binary_codes - Matrix where each row is a binary codeword
    
    if nargin < 5
        max_codes = Inf;  % Default to finding all possible codes
    end
    
    % Initialize empty code set
    binary_codes = [];
    
    % Generate all possible codewords with weight w
    % For efficiency, we only generate one representative from each cyclic equivalence class
    candidates = generate_cyclic_representatives(n, w);
    
    fprintf('Generated %d candidate codewords\n', size(candidates, 1));
    
    % For each candidate codeword
    for i = 1:size(candidates, 1)
        % Check if we've reached the maximum number of codes
        if size(binary_codes, 1) >= max_codes
            break;
        end
        
        candidate = candidates(i, :);
        
        % Check autocorrelation property
        if check_autocorrelation(candidate, lambda_a)
            
            % Check cross-correlation property with existing codes
            if check_cross_correlation(candidate, binary_codes, lambda_c)
                
                % Add candidate to code set
                binary_codes = [binary_codes; candidate];
                fprintf('Found code %d of %d\n', size(binary_codes, 1), max_codes);
            end
        end
    end
    
    fprintf('Found %d codes in total\n', size(binary_codes, 1));
end

function candidates = generate_cyclic_representatives(n, w)
    % Generate one representative from each cyclic equivalence class of weight-w binary n-tuples
    
    % Initialize to store unique canonical representatives
    canonicals = {};
    
    % Get all possible indices where 1's can be placed
    indices_combos = nchoosek(1:n, w);
    
    fprintf('Examining %d possible index combinations...\n', size(indices_combos, 1));
    
    % For each combination of indices
    for i = 1:size(indices_combos, 1)
        indices = indices_combos(i, :);
        
        % Convert to 0-indexed for easier modulo arithmetic
        indices_0based = indices - 1;
        
        % Compute the canonical form (lexicographically smallest among all cyclic shifts)
        canonical = compute_canonical_form(indices_0based, n);
        
        % Convert back to 1-indexed
        canonical = canonical + 1;
        
        % Check if this canonical form is already in our list
        is_new = true;
        for j = 1:length(canonicals)
            if isequal(canonical, canonicals{j})
                is_new = false;
                break;
            end
        end
        
        % If this is a new canonical form, add to our list
        if is_new
            canonicals{end+1} = canonical;
        end
    end
    
    fprintf('Found %d unique cyclic representatives\n', length(canonicals));
    
    % Convert canonicals to binary codewords
    candidates = zeros(length(canonicals), n);
    for i = 1:length(canonicals)
        candidates(i, canonicals{i}) = 1;
    end
end

function canonical = compute_canonical_form(indices, n)
    % Compute the canonical form of a set of indices modulo n
    % (lexicographically smallest among all cyclic shifts)
    
    % Initialize with the original indices
    canonical = sort(indices);
    
    % Check all cyclic shifts
    for shift = 1:n-1
        % Compute shifted indices
        shifted = mod(indices + shift, n);
        shifted = sort(shifted);
        
        % If this shift is lexicographically smaller, update canonical
        if lexicographically_smaller(shifted, canonical)
            canonical = shifted;
        end
    end
end

function result = lexicographically_smaller(a, b)
    % Check if array a is lexicographically smaller than array b
    
    for i = 1:min(length(a), length(b))
        if a(i) < b(i)
            result = true;
            return;
        elseif a(i) > b(i)
            result = false;
            return;
        end
    end
    
    % If we get here, the arrays are equal up to the length of the shorter one
    result = length(a) < length(b);
end

function result = check_autocorrelation(code, lambda_a)
    % Check if code satisfies autocorrelation constraint
    
    n = length(code);
    result = true;
    
    % Check for each non-zero shift
    for shift = 1:n-1
        % Calculate autocorrelation for this shift
        shifted_code = circshift(code, shift);
        correlation = sum(code & shifted_code);
        
        % If correlation exceeds lambda_a, return false
        if correlation > lambda_a
            result = false;
            return;
        end
    end
end

function result = check_cross_correlation(candidate, codes, lambda_c)
    % Check if candidate satisfies cross-correlation constraint with existing codes
    
    n = length(candidate);
    result = true;
    
    % For each existing code
    for i = 1:size(codes, 1)
        code = codes(i, :);
        
        % Check for all possible shifts
        for shift = 0:n-1
            % Calculate cross-correlation for this shift
            shifted_candidate = circshift(candidate, shift);
            correlation = sum(code & shifted_candidate);
            
            % If correlation exceeds lambda_c, return false
            if correlation > lambda_c
                result = false;
                return;
            end
        end
    end
end


function visualize_ooc(codes, lambda_a, lambda_c)
    % VISUALIZE_OOC Visualize the OOC codes and their correlation properties
    %
    % INPUTS:
    %   codes - Matrix where each row is a binary codeword
    %   lambda_a - autocorrelation constraint
    %   lambda_c - cross-correlation constraint
    
    n = size(codes, 2);
    num_codes = size(codes, 1);
    
    % Display the codes
    figure;
    subplot(2, 2, 1);
    imagesc(codes);
    title('OOC Codewords');
    xlabel('Code bit position');
    ylabel('Code index');
    colormap(gray);
    
    % Compute and display autocorrelation for first code
    if num_codes > 0
        code = codes(1, :);
        auto_corr = zeros(1, n);
        for shift = 0:n-1
            shifted_code = circshift(code, shift);
            auto_corr(shift+1) = sum(code & shifted_code);
        end
        
        subplot(2, 2, 2);
        stem(0:n-1, auto_corr, 'filled');
        title('Autocorrelation of first code');
        xlabel('Shift');
        ylabel('Correlation');
        yline(lambda_a, 'r--', ['λ_a = ' num2str(lambda_a)]);
        ylim([0 max(auto_corr)+1]);
    end
    
    % Compute and display cross-correlation between first two codes
    if num_codes > 1
        code1 = codes(1, :);
        code2 = codes(2, :);
        cross_corr = zeros(1, n);
        for shift = 0:n-1
            shifted_code2 = circshift(code2, shift);
            cross_corr(shift+1) = sum(code1 & shifted_code2);
        end
        
        subplot(2, 2, 3);
        stem(0:n-1, cross_corr, 'filled');
        title('Cross-correlation between first two codes');
        xlabel('Shift');
        ylabel('Correlation');
        yline(lambda_c, 'r--', ['λ_c = ' num2str(lambda_c)]);
        ylim([0 max(cross_corr)+1]);
    end
    
    % Display a summary of the OOC
    subplot(2, 2, 4);
    text(0.1, 0.8, ['Number of codes: ' num2str(num_codes)], 'FontSize', 12);
    text(0.1, 0.6, ['Code length (n): ' num2str(n)], 'FontSize', 12);
    text(0.1, 0.4, ['Code weight (w): ' num2str(sum(codes(1,:)))], 'FontSize', 12);
    text(0.1, 0.2, ['λ_a = ' num2str(lambda_a) ', λ_c = ' num2str(lambda_c)], 'FontSize', 12);
    axis off;
end
