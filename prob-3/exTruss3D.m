clearvars;
close all
filename='MN-P3-ExFinal-2Q-22-23.tex'; % nom del fitxer latex que generem
fout=fopen(filename,'w','n','UTF-8');  % obrir el fitxer de sortida (inclou accents)

fprintf(fout,'\\documentclass[12pt,a4paper]{article}\n'); % inici del fitxer Latex
fprintf(fout,'\\usepackage{graphicx}\n');
fprintf(fout,'\\usepackage{amsmath}\n');
fprintf(fout,'\\usepackage{moodle}\n'); % s'ha de compilar amb LuaLatex
fprintf(fout,'\\moodleset{ppi=200}\n'); 
fprintf(fout,'\\usepackage{babel}\n');  % per no tenir problemes amb accents
fprintf(fout,'\\usepackage{enumerate}\n');
fprintf(fout,'\\begin{document}\n');

category='MN-P3-ExParcial-2Q-22-23'; %nom de les preguntes com a categoria de Moodle     

fprintf(fout,'\\begin{quiz}{%s}\n',category); %(es podria no posar nom)

fprintf(fout,'\\begin{figure}[!h] \n');
fprintf(fout,'\\includegraphics[scale=0.25]{scaffold-new.png} \n');
fprintf(fout,'\\caption{Modulus of an scaffold made of cylindrical bars. \n');
fprintf(fout,'Bars $\\overline{AC}$, $\\overline{CG}$, $\\overline{GE}$, $\\overline{AE}$, \n');
fprintf(fout,'$\\overline{BD}$, $\\overline{DH}$, $\\overline{HF}$, $\\overline{BF}$, \n');
fprintf(fout,'$\\overline{AB}$, $\\overline{CD}$, $\\overline{EF}$, and $\\overline{GH}$ are of \n');
fprintf(fout,'length $L_{1} = \\ell$ mm, bars $\\overline{AG}$, and $\\overline{BH}$ \n');
fprintf(fout,'are of length $L_{2} = \\ell\\sqrt{2}$ mm, and bar $\\overline{BG}$ \n');
fprintf(fout,'is of length $L_{3} = \\ell\\sqrt{3}$ mm. }\n');
fprintf(fout,'\\end{figure} \n\n');

%The exercise has only one parameter, the value of a1 in the A in the
%initial condition g(x) = 1 + A*cos(pi*x/2)
%

% Fixed data
F = 2000.0;                 % Force modulus of the applied loads (kN)
L = 3000.0;                 % Length of Bars (mm)
Y = 212.0;                  % Young's modulus (GP, 1 GP = 1 kN/mm^2)

% Variable data
diameter = [35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0]; % (in mm), bars' 
                                                       % transversal
                                                       % section's
                                                       % dimater

nPregTotal=length(diameter); %numero de preguntes diferents que es vulguin generar.
numFakeSols = 4;

for npreg=1:nPregTotal
    % Compute the answers
    [K, Km, displ, reactForces, finalLengthOfBars] = ...
        scaffold(diameter(npreg),'n');
    
    % Specific answers & hints
    traceKm = trace(Km);                        % Answer question (a)
    K22 = K(2,2);                               % Hint question (a)

    maxRY = max(abs(reactForces(:,2)));         % Answer question (b)
    RXNod1 = reactForces(1,1);                  % Hint question (b)

    maxLengthOfBars = max(finalLengthOfBars);   % Answer question (c)
    finalLengthOfBarEM  = finalLengthOfBars(2); % hint question (c): length 
                                                % of bar EM (elem. 2 in our 
                                                % mesh)

    % Compute the fake (randomized) answers
    fake_traceKm = randomizedAnswers(traceKm,numFakeSols);
    fake_maxRY = randomizedAnswers(maxRY,numFakeSols);
    fake_maxLengthOfBars = randomizedAnswers(maxLengthOfBars,numFakeSols);

    %% Comencem a escriure la pregunta
    fprintf(fout,'\\begin{cloze}{Question %d} \n ', npreg);

    perm = ['The sketch in the Figure is part of a \n',...
         'scaffold made of cylindrical bars: $12$ bars \n',...
         'of length $L_{1} = $ %.1f mm, 2 bars of length \n',...
         '$L_{2} = L_{1}\\sqrt{2}$ mm, and $1$ bar of length \n',...
         '$L_{3} = L_{1}\\sqrt{3}$ mm, \\emph{all} them having \n',...
         '$\\phi = $ %.1f mm transversal section diameter (\\textbf{not} radius!), and \n',...
         '$Y = $ %.1f GPa Young modulus. \n\n',...
         '\\indent On the one hand we know that points $A,\\dots,F$ \n',...
         'are clamped so they can''t move in any of the three \n',...
         'spatial directions and, on the other hand, that at each of \n',...
         'the junctions $G$ and $H$ there is acting a force of modulus \n',...
         '$\\|\\vec{F}\\| = $ %.1f kN that are applied in the \n',...
         'direction of bar $\\overline{BG}$, pointing outwards the \n',...
         'structure. \n\n',...
         '\\indent As in practice 2.4, use the FEM to find the \n', ...
         'displacements of the structure''s nodes, the reaction forces, \n',...
         'and the bars'' \\emph{final} length, i.e., the bars'' length \n',...
         'when the above decribed forces are applied at the specified nodes. \n\n',...
         '\\textbf{Important remark:} set up the meshing taking the axes'' \n',...
         'position as shown in the Figure, i.e., with point B placed at the \n',...
         'origin, and the bars $\\overline{AB}$, $\\overline{BD}$, and \n',...
         '$\\overline{BF}$ on axes $x$, $y$, and $z$ respectively. \n\n',...
         'Then, answer the questions below.\n\n'
           ];
  
    fprintf(fout, perm, L, diameter(npreg), Y, F);

        primerApartat = [
            '\\begin{multi}[points=3,vertical]\n',...
            '(a) (3 points) The reduced system''s matrix''s trace''s value \n',...
            'is (in kN/mm), \n',...
            '\\item[fraction=100] %.4e\n'
            ]; 

        for e = 1:numFakeSols
            primerApartat = [
                primerApartat,...
                '\\item[fraction=-25] %.4e\n'
            ];
        end

       primerApartat = [
             primerApartat,...
             '\\item Leave it empty (no penalty) \n',...
             '\\end{multi} \n'
            ];

       primerApartat = [
           primerApartat,...
            'Hint. The entry $K_{2,2}$ of the global stiffness matrix \n',...
            'is $K_{2,2} = $ %.4e kN/mm. \n\n'
            ];
 
        fprintf(fout, primerApartat, traceKm, fake_traceKm, K22);

        segonApartat = [
            '\\begin{multi}[points=3,vertical]\n',...
             '(b) (3 points) If $Y_{j}$ denote,\n',...
             'the $y$-component of the\n',...
             'reaction force at the (global) node $j$, with $j$ ranging from \n',...
             '1 to the number of global nodes, $N$; then the absolute value \n',...
             '$\\displaystyle\\max_{j=1,\\dots,N}\\left| Y_{j}\\right|$ \n',...
             '(in kN) is\n',...
             '\\item[fraction=100] %.4e\n'
            ];

        for e = 1:numFakeSols
            segonApartat = [
                segonApartat,...
                '\\item[fraction=-25] %.4e\n'
                ];
        end

        segonApartat = [
             segonApartat,...
             '\\item Leave it empty (no penalty)\n',...
             '\\end{multi}\n\n'
             ];

         segonApartat = [
            segonApartat,...
             'Hint. The $x$-component of the reaction force of node $1$ \n',...
             'is $X_{1} = $ %.4e kN. \n\n' 
             ];
         
        fprintf(fout, segonApartat, maxRY, fake_maxRY, RXNod1);

        tercerApartat = [
            '\\begin{multi}[points=4,vertical]\n',...
            '(c) (4 points) when forces are upon, the maximum length of the \n',...
            'deformed bars (in mm) is \n',....
            '\\item[fraction=100] %.4e\n'];

        for e = 1:numFakeSols
            tercerApartat = [
                tercerApartat,...
                '\\item[fraction=-25] %.4e\n'          
                ];
        end

        tercerApartat = [
             tercerApartat,...
             '\\item Leave it empty (no penalty)\n',...
             '\\end{multi}\n\n'
             ];

        tercerApartat = [
            tercerApartat,...
             'Hint. $\\overline{CG}$ bar''s \\emph{final} length is \n',...
             '$\\left|\\overline{CG}\\right| = $ %.4e mm.\n\n'
             ];

        fprintf(fout, tercerApartat, maxLengthOfBars,...
            fake_maxLengthOfBars, finalLengthOfBarEM);

%         adviceMsg=[
%             'You can use Matlab as desk calculator and to \n',...
%             'solve the linear systems as well. \n\n'
%             ];
%         
%         fprintf(fout,adviceMsg);
  
fprintf(fout,'\\end{cloze}');
fprintf(fout,'\n\n'); 

end

fprintf(fout,'\\end{quiz}\n');
fprintf(fout,'\\end{document}  \n');

fclose(fout);

%% (c) Numerical Factory 2020 (by T. Susin, M. Calle, P. Gutierrez, J Puig)
