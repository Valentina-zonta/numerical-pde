

    meshIndex = 0;  
    topology = load(sprintf('mesh%d.topol', meshIndex)); 
    coordinates = load(sprintf('mesh%d.coord', meshIndex)); 
    boundaryConditions = load(sprintf('mesh%d.bound', meshIndex));  
    neumannBoundary = load(sprintf('mesh%d.trace', meshIndex));   
    testPoints = load(sprintf('mesh%d.track', meshIndex));   
    referenceSolution = load('solRef.dat.txt');    

    numNodes = size(coordinates,1); 
    numElements = size(topology,1); 
    timeMax = 10;  
    timeStep = 0.02; 

    % FEM implementation
    matrixH = findH(topology, coordinates, numNodes, numElements);  
    matrixM = findM(topology, coordinates, numNodes, numElements);  

    heatSolution = heatSolver(matrixM, matrixH, timeStep, timeMax, numNodes, boundaryConditions);
    [laplaceSolution, numIterations, relativeResidual] = laplaceSolver(matrixH, numNodes, boundaryConditions);

    % Interpolation and error calculation
    interpolant = scatteredInterpolant(referenceSolution(:,1), referenceSolution(:,2), referenceSolution(:,3));
    interpolatedSolution = interpolant(coordinates(:,1), coordinates(:,2));
    epsilon = laplaceError(laplaceSolution, topology, coordinates, interpolatedSolution, numNodes, numElements);

    % Mesh parameter calculation
    meshParam = meshParameter(topology, coordinates, numElements);

    % Space-time error calculation
    referenceValues = [0.2434390, 0.0772287, 0.0183037; 0.6046775, 0.2751718, 0.0928241; 0.7454968, 0.4328630, 0.1716526; 0.7751273, 0.4805514, 0.2008722];
    spaceTimeErrors = spaceTimeError(heatSolution, testPoints, timeMax, timeStep, referenceValues);
figure(1);
% graph of the steady-state solution
graph = trisurf(topology, coordinates(:,1), coordinates(:,2), laplaceSolution,...
    'EdgeColor', 'none', 'FaceColor', 'interp');
xlabel('x axis');
ylabel('y axis');
zlabel('u');
saveas(gcf, sprintf('steady_state_solution_mesh%d.png', meshIndex));

figure(2);
% graph of the solution along Neumann outer boundary at four different times
neumannSolution = zeros(size(neumannBoundary,1),4);
for i = 1:4
    tempSolution = heatSolution(:,(i*timeMax/(4*timeStep))+1);
    neumannSolution(:,i) = tempSolution(neumannBoundary(:,1));
    plot(neumannBoundary(:,2), neumannSolution(:,i), 'LineWidth', 2);
    hold on
end
legend('t = 2.5','t = 5','t = 7.5','t = 10');
ylabel('u');
xlabel('arc length');
hold off
saveas(gcf, sprintf('neumann_boundary_solution_mesh%d.png', meshIndex));

figure(3);
% graph of the solution as function of only time at the track points
trackSolution = zeros((timeMax/timeStep)+1,3);
for i = 1:3
    trackSolution(:,i) = heatSolution(testPoints(i),:);
    plot(0:timeStep:timeMax, trackSolution(:,i), 'LineWidth', 2);
    hold on
end
legend('P1','P2','P3');
ylabel('u');
xlabel('t');
hold off

figure(4);
% convergence of PCG for steady-state solution
semilogy(0:1:numIterations, relativeResidual,'-o');
xlabel('Iteration number');
ylabel('Relative residual norm');

% table 1 space-time values 
referenceValues = [0.2434390, 0.0772287, 0.0183037;...
    0.6046775, 0.2751718, 0.0928241;...
    0.7454968, 0.4328630, 0.1716526;...
    0.7751273, 0.4805514, 0.2008722];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [b,c] = computeCoeff(topology, coordinates, i, k)
    if k == 1
        b = coordinates(topology(i,2),2) - coordinates(topology(i,3),2);
        c = coordinates(topology(i,3),1) - coordinates(topology(i,2),1); 
    end

    if k == 2
        b = coordinates(topology(i,3),2) - coordinates(topology(i,1),2);
        c = coordinates(topology(i,1),1) - coordinates(topology(i,3),1); 
    end

    if k == 3
        b = coordinates(topology(i,1),2) - coordinates(topology(i,2),2);
        c = coordinates(topology(i,2),1) - coordinates(topology(i,1),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function delta = computeDelta(topology, coordinates, i)
    delta = 0.5 * det([1, coordinates(topology(i,1),1), coordinates(topology(i,1),2);...
        1, coordinates(topology(i,2),1), coordinates(topology(i,2),2);...
        1, coordinates(topology(i,3),1), coordinates(topology(i,3),2)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function matrixH = findH(topology, coordinates, numNodes, numElements)
    matrixH = zeros(numNodes, numNodes);

    for i = 1:numElements
        [b1, c1] = computeCoeff(topology, coordinates, i, 1);
        [b2, c2] = computeCoeff(topology, coordinates, i, 2);
        [b3, c3] = computeCoeff(topology, coordinates, i, 3);
        delta = computeDelta(topology, coordinates, i);

        Hloc = (1/(4*delta)) * [b1*b1+c1*c1, b1*b2+c1*c2, b1*b3+c1*c3;...
                                b1*b2+c1*c2, b2*b2+c2*c2, b2*b3+c2*c3;...
                                b1*b3+c1*c3, b2*b3+c2*c3, b3*b3+c3*c3];

        matrixH(topology(i,1), topology(i,1)) = matrixH(topology(i,1), topology(i,1)) + Hloc(1,1);
        matrixH(topology(i,1), topology(i,2)) = matrixH(topology(i,1), topology(i,2)) + Hloc(1,2);
        matrixH(topology(i,1), topology(i,3)) = matrixH(topology(i,1), topology(i,3)) + Hloc(1,3);
        matrixH(topology(i,2), topology(i,1)) = matrixH(topology(i,2), topology(i,1)) + Hloc(2,1);
        matrixH(topology(i,2), topology(i,2)) = matrixH(topology(i,2), topology(i,2)) + Hloc(2,2);
        matrixH(topology(i,2), topology(i,3)) = matrixH(topology(i,2), topology(i,3)) + Hloc(2,3);
        matrixH(topology(i,3), topology(i,1)) = matrixH(topology(i,3), topology(i,1)) + Hloc(3,1);
        matrixH(topology(i,3), topology(i,2)) = matrixH(topology(i,3), topology(i,2)) + Hloc(3,2);
        matrixH(topology(i,3), topology(i,3)) = matrixH(topology(i,3), topology(i,3)) + Hloc(3,3);
    end
    matrixH = sparse(matrixH);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function matrixM = findM(topology, coordinates, numNodes, numElements)
    matrixM = zeros(numNodes, numNodes);

    for i = 1:numElements
        delta = computeDelta(topology, coordinates, i);

        Mloc = (delta/12) * [2, 1, 1; 1, 2, 1; 1, 1, 2];

        matrixM(topology(i,1), topology(i,1)) = matrixM(topology(i,1), topology(i,1)) + Mloc(1,1);
        matrixM(topology(i,1), topology(i,2)) = matrixM(topology(i,1), topology(i,2)) + Mloc(1,2);
        matrixM(topology(i,1), topology(i,3)) = matrixM(topology(i,1), topology(i,3)) + Mloc(1,3);
        matrixM(topology(i,2), topology(i,1)) = matrixM(topology(i,2), topology(i,1)) + Mloc(2,1);
        matrixM(topology(i,2), topology(i,2)) = matrixM(topology(i,2), topology(i,2)) + Mloc(2,2);
        matrixM(topology(i,2), topology(i,3)) = matrixM(topology(i,2), topology(i,3)) + Mloc(2,3);
        matrixM(topology(i,3), topology(i,1)) = matrixM(topology(i,3), topology(i,1)) + Mloc(3,1);
        matrixM(topology(i,3), topology(i,2)) = matrixM(topology(i,3), topology(i,2)) + Mloc(3,2);
        matrixM(topology(i,3), topology(i,3)) = matrixM(topology(i,3), topology(i,3)) + Mloc(3,3);
    end
    matrixM = sparse(matrixM);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function solution = heatSolver(matrixM, matrixH, timeStep, timeMax, numNodes, boundaryConditions)
    K1 = matrixM / timeStep + matrixH / 2; 
    K2 = matrixM / timeStep - matrixH / 2;  

    A = K1;
    % imposition of Dirichlet b.c. on A
    for i = 1:size(boundaryConditions,1)
        A(:, boundaryConditions(i,1)) = 0; 
        A(boundaryConditions(i,1), :) = 0; 
        A(boundaryConditions(i,1), boundaryConditions(i,1)) = 1;
    end

    L = ichol(A, struct('type','ict','droptol',1e-3));

    tol = 1e-8;
    maxit = 1000;

    rboundary = boundaryConditions(:,2);

    solution = zeros(numNodes,(timeMax/timeStep)+1); 

    for i = 1:(timeMax/timeStep)+1
        if i > 1
            b = K2 * solution(:,i-1);
        else
            b = zeros(numNodes,1);
        end

        if i < (timeMax/(2*timeStep))+2
            for k = 1:size(boundaryConditions,1) 
                if boundaryConditions(k,2) ~= 0
                    rboundary(k) = (i-1) * 0.02 / 5;
                end
            end
        end

        % imposition of Dirichlet b.c. on b
        for k = 1:size(boundaryConditions,1)
            b = b - (K1(:,boundaryConditions(k,1)) .* rboundary(k)); 
        end

        for k = 1:size(boundaryConditions,1)
            b(boundaryConditions(k,1)) = rboundary(k); 
        end 

        % solve linear system through pcg matlab function
        [solution(:,i),~] = pcg(A, b, tol, maxit, L, L');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function epsilon = laplaceError(laplaceSolution, topology, coordinates, interpolatedSolution, numNodes, numElements)
    epsilon = zeros(numNodes,1);

    for j = 1:numElements
        delta = computeDelta(topology, coordinates, j);
        for i = 1:3
            epsilon(topology(j,i)) = epsilon(topology(j,i)) + (((interpolatedSolution(topology(j,i))) -...
                (laplaceSolution(topology(j,i))))^2) * (delta/3);
        end
    end

    epsilon = sum(epsilon);
    epsilon = sqrt(epsilon);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [solution, iter, relresvec] = laplaceSolver(matrixH, numNodes, boundaryConditions)
    A = matrixH;
    % imposition of Dirichlet b.c. on A
    for i = 1:size(boundaryConditions,1)
        A(:, boundaryConditions(i,1)) = 0; 
        A(boundaryConditions(i,1), :) = 0; 
        A(boundaryConditions(i,1), boundaryConditions(i,1)) = 1;
    end

    L = ichol(A, struct('type','ict','droptol',1e-3));

    b = zeros(numNodes,1);

    % imposition of Dirichlet b.c. on b
    for k = 1:size(boundaryConditions,1)
        b = b - (matrixH(:,boundaryConditions(k,1)) .* boundaryConditions(k,2)); 
    end

    for k = 1:size(boundaryConditions,1)
        b(boundaryConditions(k,1)) = boundaryConditions(k,2); 
    end 

    tol = 1e-8;
    maxit = 1000;

    % solve linear system through pcg matlab function
    [solution,~,~, iter, resvec] = pcg(A, b, tol, maxit, L, L');

    normb = norm(b,2);
    relresvec = resvec ./ normb;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maxSideLength = meshParameter(topology, coordinates, numElements)
    allSides = zeros(3*numElements,1);

    for i = 1:numElements
        allSides((i-1)*3+1) = norm(coordinates(topology(i,1))-coordinates(topology(i,2)),2);
        allSides((i-1)*3+2) = norm(coordinates(topology(i,1))-coordinates(topology(i,3)),2);
        allSides((i-1)*3+3) = norm(coordinates(topology(i,2))-coordinates(topology(i,3)),2);
    end

    maxSideLength = max(allSides);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spaceTimeErrors = spaceTimeError(heatSolution, testPoints, timeMax, timeStep, referenceValues)
    % Extract space-time values for the track points
    trackSolution = zeros((timeMax/timeStep)+1, 3);
    for i = 1:3
        trackSolution(:, i) = heatSolution(testPoints(i), :);
    end

    % Compute space-time values at four specific times
    computedValues = [trackSolution((timeMax/(4*timeStep))+1, 1), trackSolution((timeMax/(4*timeStep))+1, 2), trackSolution((timeMax/(4*timeStep))+1, 3);...
                      trackSolution((2*timeMax/(4*timeStep))+1, 1), trackSolution((2*timeMax/(4*timeStep))+1, 2), trackSolution((2*timeMax/(4*timeStep))+1, 3);...
                      trackSolution((3*timeMax/(4*timeStep))+1, 1), trackSolution((3*timeMax/(4*timeStep))+1, 2), trackSolution((3*timeMax/(4*timeStep))+1, 3);...
                      trackSolution((4*timeMax/(4*timeStep))+1, 1), trackSolution((4*timeMax/(4*timeStep))+1, 2), trackSolution((4*timeMax/(4*timeStep))+1, 3)];

    % Calculate space-time errors
    spaceTimeErrors = abs(computedValues - referenceValues);
end