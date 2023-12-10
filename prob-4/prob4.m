clearvars
close all
clc
fprintf('Problem 4\n')

kc = 0.05; % W/cm/C
f = 0.033; % W/cm**2
betaInt = 3.0e-3; TInfInt = 60.0;
betaOut = 5.0e-4; TInfOut = 20.0;
hint_nod_b = 117;
TcBarycenter = 67.17;
TcBarycenterHint = 67.12;

eval('curved_strip')

numNodes= size(nodes,1);
numElem= size(elem,1);

numbering= 0; %= 1 shows nodes and element numbering
plotElementsOld(nodes,elem,numbering)
hold on;

% Find Boundary points
indUpper = find(nodes(:,2) > 1.99);
indLeft = find(nodes(:,1) < 0.01);

indAboveSL = find(2.9*nodes(:,1) - 4.9*nodes(:,2)-4.41 > 0);
indNodBd = boundaryNodes(nodes, elem);
indIntNods = setdiff(indNodBd,[indAboveSL;indUpper;indLeft]);
indOutNods = setdiff(indNodBd,[indUpper;indLeft;indIntNods]);

hold on
plot(nodes(indIntNods,1),nodes(indIntNods,2),'ok','lineWidth',1,'markerFaceColor',...
    'red','markerSize',5)
plot(nodes(indOutNods,1),nodes(indOutNods,2),'ok','lineWidth',1,'markerFaceColor',...
    'blue','markerSize',5)
hold off

% Answer to part (a)
%
% Fancy output with fprintf: in exams DO NOT waste tour time with this!
%
fprintf(['(a) Sum of the x-coordinates of the pointd belonging to the top\n', ...
         '    curved boundary: %.4e\n'],sum(nodes(indIntNods,1)))
fprintf('    Hint. The sum of the y-coordinates is: %.4e\n',sum(nodes(indIntNods,2)))

%Define the coefficients vector of the model equation
a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
%f=0;
coeff=[a11,a12,a21,a22,a00,f];

%Compute the global stiff matrix
K=zeros(numNodes);    %global stiff matrix
F=zeros(numNodes,1);  %global internal forces vector
Q=zeros(numNodes,1);  %global secondary variables vector

for e = 1:numElem
    [Ke, Fe] = bilinearQuadElement(coeff,nodes,elem,e);
    rows= [elem(e,1); elem(e,2); elem(e,3); elem(e,4)];
    cols= rows;
    K(rows,cols)= K(rows,cols)+Ke;
    if (coeff(6) ~= 0)
        F(rows)= F(rows) + Fe;
    end
end %end for elements
Kini=K; %We save a copy of the initial K and F arrays
Fini=F; %for the post-process step 

%Booundary Conditions
%fixedNodes= [indB', indC'];                %fixed Nodes (global numbering)
%freeNodes= setdiff(1:numNodes,fixedNodes); %free Nodes (global numbering)

%------------- Convetion BC
[K,Q]=applyConvQuad(indIntNods,betaInt,TInfInt,K,Q,nodes,elem);
[K,Q]=applyConvQuad(indOutNods,betaOut,TInfOut,K,Q,nodes,elem);

Fm = Q + F;

%------------- Compute nodal temperature
u = K\Fm;

%------------- Post process: plot the solution
colorScale = 'jet';
titol = 'Temperature distribution';
plotContourSolution(nodes,elem,u,titol,colorScale);

% Answer to part (b)
%
% Fancy output with fprintf: in exams DO NOT waste tour time with this!
%
% fprintf('(b) The highest temperature in the domain is, T = %.5e\n', max(u))
% fprintf('Hint. Node %d has temperature T = %.5e\n',...
%     hint_nod_b, u(hint_nod_b))
fprintf('(b) The highest temperature in the domain is, T = %.4f%cC\n', ...
    max(u), char(176))
fprintf('    Hint. Node %d has temperature T = %.4f%cC\n',...
    hint_nod_b, u(hint_nod_b),char(176))

tempBarycenter = zeros(numElem,1);
for e = 1:numElem
    tempBarycenter(e) = 0.25*sum(u(elem(e,:)));
end

% Answer to part (c)
%
% Fancy output with fprintf: in exams DO NOT waste tour time with this!
%
fprintf('(c) Number of elements at temp. at barycenter > %.2f%cC: %d\n',...
    TcBarycenter, char(176), length(find(tempBarycenter > TcBarycenter)))
fprintf('    Hint. Number of elements that have a temperature < %.2f%cC: %d\n',...
    TcBarycenterHint, char(176), length(find(tempBarycenter < TcBarycenterHint)))