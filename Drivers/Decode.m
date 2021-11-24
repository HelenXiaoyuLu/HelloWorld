    query = ["I", "L", "V", "P", "M", "C", "A", "G", "P", "T", "S",...
            "Y", "W", "Q", "N", "H", "E", "D", "K", "R"];
    primers = [];
    T_degeff = getEfficientCodon();
    [fullcover, T_nonredun] = findfullcover(query, T_degeff);
    while ~any(fullcover)
        T_nonredun = sortrows(T_nonredun, 'Score', 'descend');
        primercount = 1;
        opt = T_nonredun.DegCodon(T_nonredun.Possibility == max(T_nonredun.Possibility));
        primers = [primers, opt(1)];
        encodedsubaa = T_nonredun.EncodeAA{opt(1)};
        query = setdiff(query, encodedsubaa);
        [fullcover, T_nonredun] = findfullcover(query, T_nonredun);
    end
    finaldegcodon = T_nonredun.DegCodon(fullcover);
    primers = [primers, finaldegcodon(1)];
function [fullcover, T_nonredun] = findfullcover(q, T)
    AA = ["I", "L", "V", "P", "M", "C", "A", "G", "P", "T", "S",...
        "Y", "W", "Q", "N", "H", "E", "D", "K", "R", "*"];
    qcompen = setdiff(AA, q);
    isredun = @(s) any(contains(s, qcompen));
    redunidx = cellfun(isredun, T.EncodeAA);
    T_nonredun = T(~redunidx, :);
    containsall = @(s) all(ismember(q, s)) && all(ismember(s, q));
    fullcover = cellfun(containsall, T_nonredun.EncodeAA);
end
function T = getDegDict()
    T = table();
    T.("UB Code") = ["B"; "D"; "H"; "K"; "M"; "N"; "R"; "S"; "V"; "W"; "Y"; "A"; "T"; "G"; "C"];
    T.("Description") = ["GTC"; "GAT"; "ATC"; "GT"; "AC";...
        "ACGT"; "AG"; "GC"; "ACG"; "AT"; "CT"; "A"; "T"; "G"; "C"];
    T.("Possibility") = (strlength(T.Description) + 1)./2;
    T.Properties.RowNames = T.("UB Code");
end
function T = getEncodingDict()
    AA = ["I", "L", "V", "F", "M", "C", "A", "G", "P", "T", "S",...
        "Y", "W", "Q", "N", "H", "E", "D", "K", "R", "*"];
    Codon = [{{"ATT", "ATC", "ATA"}}, ...
        {{"CTT", "CTC", "CTA", "CTG", "TTA", "TTG"}}, ...
        {{"GTT", "GTC", "GTA", "GTG"}}, ...
        {{"TTT", "TTC"}}, ...
        {{"ATG"}}, ...
        {{"TGT", "TGC"}}, ...
        {{"GCT", "GCC", "GCA", "GCG"}}, ...
        {{"GGT", "GGC", "GGA", "GGG"}}, ...
        {{"CCT", "CCC", "CCA", "CCG"}}, ...
        {{"ACT", "ACC", "ACA", "ACG"}}, ...
        {{"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}}, ...
        {{"TAT", "TAC"}}, ...
        {{"TGG"}}, ...
        {{"CAA", "CAG"}}, ...
        {{"AAT", "AAC"}}, ...
        {{"CAT", "CAC"}}, ...
        {{"GAA", "GAG"}}, ...
        {{"GAT", "GAC"}}, ...
        {{"AAA", "AAG"}}, ...
        {{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}}, ...
        {{"TAA", "TAG", "TGA"}}];
    T = dict(AA, Codon);
end
function score = calcscore(aa)
    T = getEncodingDict();
    score = 0;
    for i = 1:numel(aa)
        score = score + 1/numel(T{aa(i)}{:});
    end
end
function T_degeff = getEfficientCodon()
    T = getDegDict();
    degcodon = {};
    encodeaas = {};
    effecient = [];
    degcodoncounter = 1;
    for i = 1:height(T)
        for j = 1:height(T)
            for k = 1:height(T)
                degcodon{degcodoncounter} = T.("UB Code")(i) + T.("UB Code")(j) + T.("UB Code")(k);
                undegcodon = [T.Description(i), T.Description(j), T.Description(k)];
                c1 = char(undegcodon(1));
                c2 = char(undegcodon(2));
                c3 = char(undegcodon(3));
                ncodon = [length(c1), length(c2), length(c3)];
                encodeaa = strings(1, ncodon(1)*ncodon(2)*ncodon(3));
                aacounter = 1;
                for u = 1:ncodon(1)
                    for v = 1:ncodon(2)
                        for w = 1:ncodon(3)
                            singcodon = [c1(u), c2(v), c3(w)];   
                            encodeaa(aacounter) = nt2aa(singcodon, 'ALTERNATIVESTARTCODONS', false);
                            aacounter = aacounter + 1;
                        end
                    end
                end
                encodeaas{degcodoncounter} = encodeaa;
                if numel(unique(encodeaa)) < numel(encodeaa) % replicated encoding
                    effecient(degcodoncounter) = 0;
                else
                    effecient(degcodoncounter) = 1;
                end
                degcodoncounter = degcodoncounter + 1;
            end
        end
    end

    T_deg = table(cat(1,degcodon{:}), encodeaas', effecient', ...
        'VariableNames', ["DegCodon", "EncodeAA", "Efficient"]);
    T_degeff = T_deg(T_deg.Efficient == 1, :);
    T_degeff.Properties.RowNames = T_degeff.DegCodon;
    T_degeff.("Possibility") = cellfun(@length, T_degeff.EncodeAA);
    T_degeff.("Score") = cellfun(@calcscore, T_degeff.EncodeAA);
end
% function n = numericRepresent(degcodon)
%     A = 1;
%     T = 2;
%     G = 4;
%     C = 8;
%     B = G+T+C;
%     D = G+A+T;
%     H = A+T+C;
%     K = G+T;
%     M = A+C;
%     N = A+T+G+C;
%     R = A+G;
%     S = G+C;
%     V = A+C+G;
%     W = A+T;
%     Y = C+T;
% end