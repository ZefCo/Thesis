%  Program to convert the intron and exon sequences to .png images.
%  Points are `scatterred' and images are nomalized to sum=1.

%InputFile   = '../Data/Intron.xlsx';
%OutputDir   = '../Images/Images_Q3_3K/Intron/';
InputFile   = 'Intron.xlsx';
OutputDir   = 'Intron/';

LattSize    = 64;
QuartLength = 3;
EdgeLength  = 10;
NSide       = 4^QuartLength;
NumImages   = 3000;

% Factors and Sums for Fuzzy Distribution
Fact1  = exp(-0.5);            Fact2   = exp(-1);
SumCent = 1+4*(Fact1+Fact2);   SumEdge = 1+3*Fact1+2*Fact2;   SumVert = 1+2*Fact1+Fact2;

% Assignment of numbers to symbols - Can be Different for Anterior and
% Posterior Sub-Sequences
Post('A')=0; Post('G')=1; Post('T') = 2; Post('C')=3;
Post('a')=0; Post('g')=1; Post('t') = 2; Post('c')=3;
Antr('A')=0; Antr('G')=1; Antr('T') = 2; Antr('C')=3;
Antr('a')=0; Antr('g')=1; Antr('t') = 2; Antr('c')=3;

% Read Data
AllDat    = readtable (InputFile);
Posterior = table2array (AllDat (:, 3));
Sequences = table2array (AllDat (:, 4));
Anterior  = table2array (AllDat (:, 5));

for nseq = 1:1:NumImages %size(Sequences, 1)
    clear x y Density;

    Strn = strcat (cell2mat(Posterior(nseq)), cell2mat(Sequences(nseq)), cell2mat(Anterior(nseq)));
    Density = zeros(NSide);                 % Density Array Reset

    for np = (EdgeLength+1):1:(EdgeLength+numel(cell2mat(Sequences(nseq))))
        x(np) = 0;
        y(np) = 0;
        for nq = 1:1:QuartLength
            x(np) = x(np) + (Antr(Strn(np+nq-1)) / 4^nq);
            y(np) = y(np) + (Post(Strn(np-nq)) / 4^nq);
        end

        Lattx = NSide * x(np);
        Latty = NSide * y(np);

        if (Lattx == 0)
            if (Latty == 0)
                Density(1, 1) = Density(1, 1) + 1/SumVert;
                Density(1, 2) = Density(1, 2) + Fact1/SumVert;
                Density(2, 1) = Density(2, 1) + Fact1/SumVert;
                Density(2, 2) = Density(2, 2) + Fact2/SumVert;
            elseif (Latty < NSide-1)
                Density(1, Latty+1) = Density(1, Latty+1) + 1/SumEdge;
                Density(1, Latty)   = Density(1, Latty) + Fact1/SumEdge;
                Density(1, Latty+2) = Density(1, Latty+2) + Fact1/SumEdge;
                Density(2, Latty+1) = Density(2, Latty+1) + Fact1/SumEdge;
                Density(2, Latty)   = Density(2, Latty) + Fact2/SumEdge;
                Density(2, Latty+2) = Density(2, Latty+2) + Fact2/SumEdge;
            else
                Density(1, NSide)   = Density(1, NSide) + 1/SumVert;
                Density(1, NSide-1) = Density(1, NSide-1) + Fact1/SumVert;
                Density(2, NSide)   = Density(2, NSide) + Fact1/SumVert;
                Density(2, NSide-1) = Density(2, NSide-1) + Fact2/SumVert;
            end
        elseif (Lattx < NSide-1)
            if (Latty == 0)
                Density(Lattx+1, 1) = Density(Lattx+1, 1) + 1/SumEdge;
                Density(Lattx, 1)   = Density(Lattx, 1) + Fact1/SumEdge;
                Density(Lattx+2, 1) = Density(Lattx+2, 1) + Fact1/SumEdge;
                Density(Lattx+1, 2) = Density(Lattx+1, 2) + Fact1/SumEdge;
                Density(Lattx, 2)   = Density(Lattx, 2) + Fact2/SumEdge;
                Density(Lattx+2, 2) = Density(Lattx+2, 2) + Fact2/SumEdge;
            elseif (Latty < NSide-1)
                Density(Lattx+1, Latty+1) = Density(Lattx+1, Latty+1) + 1/SumCent;
                Density(Lattx+1, Latty)   = Density(Lattx+1, Latty) + Fact1/SumCent;
                Density(Lattx+1, Latty+2) = Density(Lattx+1, Latty+2) + Fact1/SumCent;
                Density(Lattx, Latty+1)   = Density(Lattx, Latty+1) + Fact1/SumCent;
                Density(Lattx+2, Latty+1) = Density(Lattx+2, Latty+1) + Fact1/SumCent;
                Density(Lattx+1, Latty)   = Density(Lattx+1, Latty) + Fact2/SumCent;
                Density(Lattx+1, Latty+2) = Density(Lattx+1, Latty+2) + Fact2/SumCent;
                Density(Lattx, Latty+1)   = Density(Lattx, Latty+1) + Fact2/SumCent;
                Density(Lattx+2, Latty+1) = Density(Lattx+2, Latty+1) + Fact2/SumCent;
            else
                Density(Lattx+1, NSide)   = Density(Lattx+1, NSide) + 1/SumEdge;
                Density(Lattx, NSide)     = Density(Lattx, NSide) + Fact1/SumEdge;
                Density(Lattx+2, NSide)   = Density(Lattx+2, NSide) + Fact1/SumEdge;
                Density(Lattx+1, NSide-1) = Density(Lattx+1, NSide-1) + Fact1/SumEdge;
                Density(Lattx, NSide-1)   = Density(Lattx, NSide-1) + Fact2/SumEdge;
                Density(Lattx+2, NSide-1) = Density(Lattx+2, NSide-1) + Fact2/SumEdge;
            end
        else
            if (Latty == 0)
                Density(NSide, 1)   = Density(NSide, 1) + 1/SumVert;
                Density(NSide, 2)   = Density(NSide, 2) + Fact1/SumVert;
                Density(NSide-1, 1) = Density(NSide-1, 1) + Fact1/SumVert;
                Density(NSide-1, 2) = Density(NSide-1, 2) + Fact2/SumVert;
            elseif (Latty < NSide-1)
                Density(NSide, Latty+1)   = Density(NSide, Latty+1) + 1/SumEdge;
                Density(NSide, Latty)     = Density(NSide, Latty) + Fact1/SumEdge;
                Density(NSide, Latty+2)   = Density(NSide, Latty+2) + Fact1/SumEdge;
                Density(NSide-1, Latty+1) = Density(NSide-1, Latty+1) + Fact1/SumEdge;
                Density(NSide-1, Latty)   = Density(NSide-1, Latty) + Fact2/SumEdge;
                Density(NSide-1, Latty+2) = Density(NSide-1, Latty+2) + Fact2/SumEdge;
            else
                Density(NSide, NSide)     = Density(NSide, NSide) + 1/SumVert;
                Density(NSide-1, NSide)   = Density(NSide-1, NSide) + Fact1/SumVert;
                Density(NSide, NSide-1)   = Density(NSide, NSide-1) + Fact1/SumVert;
                Density(NSide-1, NSide-1) = Density(NSide-1, NSide-1) + Fact2/SumVert;
            end
        end
    end
    Density = Density' / numel (cell2mat(Sequences(nseq)));

    % Save field of size LattSizeXLattSize
    % save (strcat (OutputDir, "Field", sprintf('%05d',nseq), ".mat"), "Density", "NSide");

    % Save PNG files
    contourf ([0:NSide-1]/NSide, [0:NSide-1]/NSide, Density, 25, 'LineStyle', 'none'), axis('square');
    colormap (gray(256));
    set(gca, 'XTick', []); set(gca, 'YTick', []);
    ax1.XAxis.Visible = 'off';  ax1.YAxis.Visible = 'off'; 
    exportgraphics (gcf, strcat(OutputDir, "Image", sprintf('%05d',nseq), ".png"), 'Resolution', 50);
end
