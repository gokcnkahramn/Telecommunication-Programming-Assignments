% Note: To open those questions, you need to use "Run Sections".
%% Q1
symbols = 1:5; %This is our alphabet.
p = [1/2 1/4 1/12 1/12 1/12]; % This is our pmf.
dictt = cell(5,2); % This is our dictionary for encoding.
dictt{1,1} = 1;
dictt{2,1} = 2;
dictt{3,1} = 3;
dictt{4,1} = 4;
dictt{5,1} = 5;
sum = 0; % This thing is for comparing part in Huffman encoding.
k = 4; % Since we make four comparisons in our code, this number is for comparison.
while k >= 1
    sum = p(k + 1) + p(k);
    if k == 4 % This is for first comparison.
        if p(k + 1) >= p(k)
                dictt{k + 1,2}(k) = 1;
                dictt{k,2}(k) = 0;
        else
                dictt{k,2}(k) = 1;
                dictt{k + 1,2}(k) = 0;
        end
    elseif k == 3 % This is for second comparison.
        if sum >= p(k)
                dictt{k + 2,2}(k) = 1;
                dictt{k + 1,2}(k) = 1;
                dictt{k,2}(k) = 0;
        else
                dictt{k + 2,2}(k) = 1;
                dictt{k,2}(k) = 1;
                dictt{k + 1,2}(k) = 0;
        end
    elseif k == 2 % This is for third comparison.
        if sum >= p(k)
                dictt{k + 3,2}(k) = 1;            
                dictt{k + 2,2}(k) = 1;
                dictt{k + 1,2}(k) = 1;
                dictt{k,2}(k) = 0;
        else
                dictt{k + 3,2}(k) = 1;            
                dictt{k + 2,2}(k) = 1;
                dictt{k + 1,2}(k) = 0;
                dictt{k,2}(k) = 1;
        end
    elseif k == 1 % This is for last comparison.
        if sum >= p(k)
            dictt{k + 4,2}(k) = 1;            
            dictt{k + 3,2}(k) = 1;            
            dictt{k + 2,2}(k) = 1;
            dictt{k + 1,2}(k) = 1;
            dictt{k,2}(k) = 0;
        else
            dictt{k + 4,2}(k) = 1;            
            dictt{k + 3,2}(k) = 1;            
            dictt{k + 2,2}(k) = 1;
            dictt{k + 1,2}(k) = 0;
            dictt{k,2}(k) = 1;
        end
    end
k = k - 1;
end


a = 1;
entropy = 0;
while a < 6 % This is for calculating entropy.
    entropy = entropy + length(dictt{a,2})*p(a);
    a = a + 1;
end
inputSig = randsrc(10000,1,[symbols;p]); % This is for our random symbols
fixedentropy = fix(entropy);
% This converts original signal to binary and it shows its lenght.
binarySig = int2bit(inputSig,fixedentropy);
seqLen = numel(binarySig); 
% This converts Hufmann encoded to binary, lenght for binary symbols.
aa = 1;
encodedLen = 0;
while aa < 10000
    dictelement = dictt{inputSig(aa),2};
    encodedLen = encodedLen + length(dictelement);
aa = aa + 1;    
end

%% Q2

% As I wrote report, I used "huffmanenco" algorithm.
% Because when I started to write same algorithm with question 1,
% There will be 24 different branches which is very long.
% For decreasing my work, I used "huffmanenco".
% ! Otherwise, if we do it with the previous question, it will give the same result. !
% https://www.mathworks.com/help/comm/ref/huffmanenco.html
symbols = [11 12 13 14 15 21 22 23 24 25 31 32 33 34 35 41 42 43 44 45 51 52 53 54 55]; % This is our alphabet.
p = [1/4 1/8 1/24 1/24 1/24 1/8 1/16 1/48 1/48 1/48 1/24 1/48 1/144 1/144 1/144 1/24 1/48 1/144 1/144 1/144 1/24 1/48 1/144 1/144 1/144]; % This is its probability mass function.
dict = huffmandict(symbols,p); % This is our huffman dictionary.
inputSig = randsrc(10000,1,[symbols;p]); % This generates random symbols.
encd = huffmanenco(inputSig,dict); % This makes huffman encoding
% To calculate total entropy
entropy = 0;
for i = 1:length(symbols)
    addedentropy = p(i)*length(dict{i, 2});
    entropy = entropy + addedentropy;
end
fixedentropy = fix(entropy);
% Making our huffman encoding - decoding true
sig = huffmandeco(encd,dict);
equality = isequal(inputSig,sig);
% This converts original signal to binary and it shows its lenght.
binarySig = int2bit(inputSig,fixedentropy);
seqLen = numel(binarySig); 
% This converts Hufmann encoded to binary, lenght for binary symbols.
binaryComp = int2bit(encd,fixedentropy);
encodedLen = numel(binaryComp);

%% Q3   

symbols = 1:5; %This is our alphabet.
p1 = [1/2 1/4 1/12 1/12 1/12]; % This is our pmf1.
p2 = [0.48 0.22 0.1 0.1 0.1]; % This is our pmf1.
dictt = cell(5,2); % This is our dictionary for encoding.
dictt{1,1} = 1;
dictt{2,1} = 2;
dictt{3,1} = 3;
dictt{4,1} = 4;
dictt{5,1} = 5;
sum = 0; % This thing is for comparing part in Huffman encoding.
k = 4; % Since we make four comparisons in our code, this number is for comparison.
while k >= 1
    sum = p2(k + 1) + p2(k);
    if k == 4 % This is for first comparison.
        if p2(k + 1) >= p2(k)
                dictt{k + 1,2}(k) = 1;
                dictt{k,2}(k) = 0;
        else
                dictt{k,2}(k) = 1;
                dictt{k + 1,2}(k) = 0;
        end
    elseif k == 3 % This is for second comparison.
        if sum >= p2(k)
                dictt{k + 2,2}(k) = 1;
                dictt{k + 1,2}(k) = 1;
                dictt{k,2}(k) = 0;
        else
                dictt{k + 2,2}(k) = 1;
                dictt{k,2}(k) = 1;
                dictt{k + 1,2}(k) = 0;
        end
    elseif k == 2 % This is for third comparison.
        if sum >= p2(k)
                dictt{k + 3,2}(k) = 1;            
                dictt{k + 2,2}(k) = 1;
                dictt{k + 1,2}(k) = 1;
                dictt{k,2}(k) = 0;
        else
                dictt{k + 3,2}(k) = 1;            
                dictt{k + 2,2}(k) = 1;
                dictt{k + 1,2}(k) = 0;
                dictt{k,2}(k) = 1;
        end
    elseif k == 1 % This is for last comparison.
        if sum >= p2(k)
            dictt{k + 4,2}(k) = 1;            
            dictt{k + 3,2}(k) = 1;            
            dictt{k + 2,2}(k) = 1;
            dictt{k + 1,2}(k) = 1;
            dictt{k,2}(k) = 0;
        else
            dictt{k + 4,2}(k) = 1;            
            dictt{k + 3,2}(k) = 1;            
            dictt{k + 2,2}(k) = 1;
            dictt{k + 1,2}(k) = 0;
            dictt{k,2}(k) = 1;
        end
    end
k = k - 1;
end


a = 1;
entropy = 0;
while a < 6 % This is for calculating entropy.
    entropy = entropy + length(dictt{a,2})*p2(a);
    a = a + 1;
end
inputSig = randsrc(10000,1,[symbols;p1]); % This is for our random symbols
fixedentropy = fix(entropy);
% This converts original signal to binary and it shows its lenght.
binarySig = int2bit(inputSig,fixedentropy);
seqLen = numel(binarySig); 
% This converts Hufmann encoded to binary, lenght for binary symbols.
aa = 1;
encodedLen = 0;
while aa < 10000
    dictelement = dictt{inputSig(aa),2};
    encodedLen = encodedLen + length(dictelement);
aa = aa + 1;    
end

%% Q4
binaries = 0:1; % This is our binaries. 
prob = [0.9 0.1]; % This is our probabilities. 
randlist = randsrc(1000000,1,[binaries;prob]); % This is our random list. If you want to check, you need to write 100.000 to 1.000.000
randlist(end+1) = "#"; % For all Lempel Ziv algorithm, we need to have hashtag for ending.
dict = ["#"; 0; 1]; % This is our dictionary.
nxtdict = length(dict) + 1; % This is our adding dictionary for checking next variables.
k = 1;
cont = ""; % This side is for controlling digits for adding new element to dictionary.
cont = cont + randlist(1);
% This side is for encoding.
encd = zeros(1,length(randlist));
incdt = 1;
while k < length(randlist)
    for m = 1:length(dict)
        if cont == dict(m)
            add = m;
            nextt = "";
            nextt = nextt + randlist(k+1);
            if nextt == '#'
                k = length(randlist)+1;
                break
            end
            cont = cont + nextt;
            k = k+1;
            m = 1;
        end
    end
    if nextt ~= '#'
        encd(incdt) = add;
        incdt = incdt + 1;
        dict(nxtdict) = cont;
        cont = nextt;
        nxtdict = nxtdict + 1;
        
    else
        encd(incdt) = add;
        incdt = incdt + 1;
    end
end
encd = nonzeros(encd) - 1;

% Since we convert all digits zero and one for dictionary, this is for encoding part.
tot = length(de2bi(length(dict)));
binaryencd = zeros(length(encd), tot);
for k = 1:length(encd)
    binaryencd(k,:) =  [zeros(1,tot-length(de2bi(encd(k)))) de2bi(encd(k), 'left-msb')];
end
% For getting our values, you can use those commands which names are:
% disp(dict), disp(randlist) and disp(encd)

%% Q5

p0 = 0.9;  % This is our probability for 0
p1 = 0.1;  % This is our probability for 1
N = 1000000;  % This is number of source output, you can change it for 100.000
totdict = 2;  % This is our initial dictionary.

% This is for generating random samples
src = rand(1,N) > p0;

% This is for Lempel-Ziv-Welch encoding part
dict = containers.Map({'0','1'},{0,1});  % This is our dictionary part.
encd = [];  % This is our encoded message
seq = '';  % This is our current check part
for i = 1:N
    nxt = num2str(src(i));  % This is for next output
    if isKey(dict, [seq nxt])  % This is for putting our variables in dictionary
        seq = [seq nxt];
    else  % This is for not part of dictionary
        encd(end+1) = dict(seq);  % Adding new variable
        dict([seq nxt]) = totdict;  % Add new entry to dictionary
        totdict = totdict + 1;  % Increase our dictionary by one
        seq = nxt;  % Start next sequence
    end
end
encd(end+1) = dict(seq);  % This is for last sequence which is "#"