function hw6_q1(initialState, T, P, lambdas)

%% pre-defined variables for testing purpose, need to comment

N = 10;
%{
P_unNormalized = rand(N + 1, N + 1);
sum_rows = sum(P_unNormalized,2);
sum_rows_extended = repmat(sum_rows, 1, N + 1);

% normlized to get transition probability matrix
P = P_unNormalized./sum_rows_extended;
lambdas = rand(N + 1, 1);

%}
disp('Transition Probability Matrix');
disp(P);

times = [];
states = [];
t = 0;
s = initialState;
states(end+1) = s; 
times(end+1) = t; 

t = exprnd(lambdas(s+1));

while t < T
    
    
    p = P(s + 1,:);
    
    u = rand;
    %disp(u);
    
    for i = 1:N+1
        if i == 1
            if u <= p(i)
                j = i;
                break;
            end
        elseif i == N + 1
            j = i;
            break;
        else
            if ((u > sum(p(1:i-1))) && (u <= sum(p(1:i))))
                j = i;
                break;
            end    
        end
    end
    % map state to 0 ~ N
    j = j - 1;
    disp(['U: ', num2str(u), '; t: ', num2str(t),'; transiting from state ', num2str(s), ' to state ', num2str(j)]);
    s = j;
    states(end+1) = s; 
    times(end+1) = t; 
    %pause(0.001);
    
    t = t + exprnd(lambdas(s+1));
end

%plot(times, states);
%title(['Markov Chain with ', num2str(N + 1), ' states']); 

fprintf('\nStates:\n');
fprintf('%d ', states);
fprintf('\nTimes:\n');
fprintf('%d ', times);  
fprintf('\n');
end
