clearvars
close all
clc

% ============================ Scaffold ===================================
fileSols = 'solutionEx3.txt';
fOut = fileSols;
fOut=fopen(fileSols,'w');
fprintf(fOut,'Problem 3\n\n');

% Data
F = 2000.0;                   % Force (kN)
L = 3000.0;                   % Length of Bars (mm)
diameterOfSectionBars = 35.0; % diameter of section bars (mm)
Y = 212.0;                    % Young's modulus (GP, 1 GP = 1 kN/mm^2)
%nodDisp = 4;                 % Node disp.


fprintf(fOut,'Length................................... L = %e (mm)\n',  L);
fprintf(fOut,'Modulus of force at nodes M, N, and O.... F = %e (kN)\n',  F);
fprintf(fOut,'Diameter of section bars............... Phi = %e (mm)\n',  diameterOfSectionBars);
fprintf(fOut,'Young modulus............................ Y = %e (GPa)\n', Y);
fprintf(fOut,'\n');

%Geometry
nodes = [
    0,0,1;  % node 1 (node A) 
    0,0,0;  % node 2 (node B)
    1,0,1;  % node 3 (node C)
    1,0,0;  % node 4 (node D)
    0,1,1;  % node 5 (node E)
    0,1,0;  % node 6 (node F)
    1,1,1;  % node 7 (node G)
    1,1,0   % node 8 (node H)
    ];

nodes = L*nodes;

elem = [
    1, 3;   % elem 1
    3, 7;   % elem 2
    7, 5;   % elem 3
    1, 5;   % elem 4
    2, 4;   % elem 5
    4, 8;   % elem 6
    6, 8;   % elem 7
    2, 6;   % elem 8
    1, 2;   % elem 9
    3, 4;   % elem 10
    7, 8;   % elem 11
    5, 6;   % elem 12
    1, 7;   % elem 13
    2, 7;   % elem 14
    2, 8;   % elem 15
    ];  

numbering=0; %1
plotElementsOld(nodes,elem,numbering)
%plotElements(nodes,elem,numbering)
    
%%
numNod=size(nodes,1);
numElem=size(elem,1);
ndim=size(nodes,2);

%Real constants
Area = pi*(diameterOfSectionBars/2)^2;
A = Area*ones(1, numElem);
E = Y*ones(1, numElem);

%Assembly
u=zeros(ndim*numNod, 1);
Q=zeros(ndim*numNod, 1);
K=zeros(ndim*numNod);

for e=1:numElem
    Ke=spatialLinkStiffMatrix(nodes,elem,e,E,A);
    rows=[ndim*elem(e,1) - 2, ndim*elem(e,1) - 1, ndim*elem(e,1),...
          ndim*elem(e,2) - 2, ndim*elem(e,2) - 1, ndim*elem(e,2)];
    cols=rows;
    K(rows,cols)=K(rows,cols)+Ke; %Assembly
end

% Boundary conditions
fixedNods = 1:6;
fixedNods = [ndim*fixedNods-2, ndim*fixedNods-1, ndim*fixedNods];
freeNods = setdiff(1:ndim*numNod, fixedNods);

% Natural B.C.
% Side Loads
vect = nodes(7,:) - nodes(2,:);
vect = vect/norm(vect);
R = F*vect;

%node 7
nod = 7;
Q(ndim*nod-2) = R(1);   %kN
Q(ndim*nod-1) = R(2);   %kN
Q(ndim*nod) = R(3);     %kN 

%node 8
nod = 8;
Q(ndim*nod-2) = R(1);   %kN
Q(ndim*nod-1) = R(2);   %kN
Q(ndim*nod) = R(3);     %kN

%EssentialBoundary Conditions
u(fixedNods) = 0.0;     % Not necessary...

%Reduced system
Qm = Q(freeNods) - K(freeNods,fixedNods)*u(fixedNods);
Km = K(freeNods,freeNods);

%Solve the reduced system
um = Km\Qm;
u(freeNods) = um;

% -------------------------------------------------------------------------
% Post process 
% -------------------------------------------------------------------------
% 1. Plot deformed structure
esc = 1;
plotDeformedTruss(nodes,elem,u,esc);

% 2. Displacements
displ = [u(1:ndim:end), u(2:ndim:end), u(3:ndim:end)];

% 3. Reaction forces
Fr = K*u - Q;
reactForces = [Fr(1:ndim:end), Fr(2:ndim:end), Fr(3:ndim:end)];

% 4. Final length of bars
finalLengthOfBars = zeros(numElem, 1);
for e = 1:numElem
    n1 = elem(e,1);
    n2 = elem(e,2);
    finalLengthOfBars(e) = ...
        norm(nodes(n2,:) + displ(n2,:) - nodes(n1,:) - displ(n1,:));
end

nodDispl = (1:numNod)';

% Print out all the displacements at each node
fprintf(fOut, '\n%21s\n', 'Displacements (in mm)');
fprintf(fOut, '%9s%11s%17s%17s\n', 'NOD.','UX','UY','UZ');
fprintf(fOut, '%7d%17.5e%17.5e%17.5e\n',[nodDispl, displ]');

% Print out all the reaction forces at each node
fprintf(fOut, '\n%23s\n', 'Reaction forces (in kN)');
fprintf(fOut, '%9s%11s%17s%17s\n', 'NOD.','RX', 'RY', 'RZ');
fprintf(fOut, '%7d%17.5e%17.5e%17.5e\n', [nodDispl, reactForces]');

% Print out final lenght of bars (after nodal forces are applied)
elemLength = (1:numElem)';
fprintf(fOut, '\n%32s\n', 'Final length of the bars (in mm)');
fprintf(fOut, '%9s%14s\n', 'ELEM.', 'FINAL LENGTH');
fprintf(fOut, '%7d%16.5e\n', [elemLength, finalLengthOfBars]');

% Solutions
fprintf(fOut, '\nSolutions:\n');
fprintf(fOut, '(a) Reduced system''s matrix''s trace, tr(Km) = %.4e\n', trace(Km));
fprintf(fOut, '    Hint. K(2,2) = %.4e\n', K(2,2));
fprintf(fOut, '(b) max |Y(i)| = %.4e\n', max(abs(reactForces(:,2))));
fprintf(fOut, '    Hint. x-component of the reaction force at node 1, RX(1) = %.4e kN\n', ...
    reactForces(1,1));
fprintf(fOut, '(c) Bars'' maximum length once the forces are applied, maxLength = %.4e mm\n', ...
    max(finalLengthOfBars));
fprintf(fOut, '    Hint. CG bar''s final length, length(CG) = %.4e mm\n', ...
    finalLengthOfBars(2));

fclose(fOut);
type(fileSols);
