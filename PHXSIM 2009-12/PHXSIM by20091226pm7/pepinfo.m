%%%2009-12-10 pepinfo.m: 

function [peptideMass, peptideCarbonNum]=pepinfo(subSeq)

AAshort  =  ['A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y'];

AAcarbonNum=[ 3,  3,  4,  5,  9,  2,  6,  6,  6,  6,  5,  4,  12, 5,  5,  6,  3,  4,  3,  5,  11, 9 ];

AAmonoMass=[71.037110
103.009190
115.026940
129.042590
147.068410
57.021460
137.058910
113.084060
128.094960
113.084060
131.040490
114.042930
255.158290
97.052760
128.058580
156.101110
87.032030
101.047680
168.964200
99.068410
186.079310
163.063330]; %above values from http://en.wikipedia.org/wiki/Proteinogenic_amino_acid


peptideMass=0;
peptideCarbonNum=0;

for i=1:size(subSeq,2)
    index=find(subSeq(i)==AAshort);
    peptideMass=peptideMass+AAmonoMass(index);
    peptideCarbonNum=peptideCarbonNum+AAcarbonNum(index);
end

peptideMass=peptideMass + (1.007825*2+15.994915); %peptide's mass is the sum of the residue masses plus the mass of water.



        

        