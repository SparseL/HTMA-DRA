function f = tournament_selection(chromosome, pool_size, tour_size, fit)

% Get the size of chromosome. The number of chromosome is not important
% while the number of elements in chromosome are important.
[pop, ~] = size(chromosome);
% Until the mating pool is filled, perform tournament selection
for i = 1 : pool_size
    % Select n individuals at random, where n = tour_size
    for j = 1 : tour_size
        % Select an individual at random
        candidate(j) = round(pop*rand(1));
        % Make sure that the array starts from one. 
        if candidate(j) == 0
            candidate(j) = 1;
        end
        if j > 1
            % Make sure that same candidate is not choosen.
            while ~isempty(find(candidate(1 : j - 1) == candidate(j)))
                candidate(j) = round(pop*rand(1));
                if candidate(j) == 0
                    candidate(j) = 1;
                end
            end
        end
    end
    % Collect information about the selected candidates.
    for j = 1 : tour_size
        c_obj(j) = fit(candidate(j));
    end
    % Find the candidate with the least rank
    min_candidate = ...
        find(c_obj == min(c_obj));
    % Add the selected individual to the mating pool
    f(i,:) = chromosome(candidate(min_candidate(1)),:);
end
