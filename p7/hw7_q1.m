m = 1000;
qArray = zeros(m+1,2);
qArray(1,:) = [0.5 0.50];

e = 0.05;
L = 20;

for i = 2:m
    qCurrent = qArray(i-1,:);
    q = qCurrent;
    d = length(q);

    p = randn(1,d);
    pCurrent = p;
    
    for j = 1:L
        q = q + e*p; % make a full step for the position
    end
    
    % To make the proposal symmetric, negate momentum at the end of trajectory
    p = -p;

    % Evaluate potential & kinetic energies to start and end the trajectory
    U = hw7_q1_evaluateU(qCurrent(1), qCurrent(2));
    K = sum(pCurrent.^2)/2;

    % proposed U
    U_tilde = hw7_q1_evaluateU(q(1), q(2));
    % proposed K
    K_tilde = sum(p.^2)/2;
    
    % accept/reject the state at the end of trajectory.
    % i.e. returning either the position at the end of the trajectory or the initial position
    
    u = rand;
    if (U_tilde == 0) || (u > exp( U-U_tilde + K-K_tilde))
        qArray(i,:) = qCurrent;
    else
        qArray(i,:) = q;
    end

end

q1Array = qArray(:,1);
q2Array = qArray(:,2);

%% Traceplots

figure(1);

plot([-1 -1],[0 -1], 'b', 'LineWidth', 2);
hold on;
plot([1 1],[0 1], 'b', 'LineWidth', 2);
plot([-1 1],[0 0], 'b', 'LineWidth', 2);
plot([-1 0],[-1 -1], 'b', 'LineWidth', 2);
plot([1 0],[1 1], 'b', 'LineWidth', 2);
plot([0 0],[-1 1], 'b', 'LineWidth', 2);

plot(q1Array, q2Array);
xlabel('x');
ylabel('y')

hold off;

figure(2);

scatter(q1Array, q2Array)
