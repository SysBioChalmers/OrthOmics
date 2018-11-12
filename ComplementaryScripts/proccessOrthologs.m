function [uniprots,genes] = proccessOrthologs
%   proccessOrthologs
%
%   Function that, given a OrthoFinder results file, generates tables with 
%   Ortho groups IDs and the respective proteins and genes IDs for each of
%   the analyzed organisms.
%
%   Ivan Domenzain.     2018-09-19
%

%%Load gene-protein databases
load('../Databases/kma_ProtDatabase.mat')
kma_kegg = kegg;
kma_swissprot = swissprot;
load('../Databases/yli_ProtDatabase.mat')
yli_kegg = kegg;
yli_swissprot = swissprot;
load('../Databases/sce_ProtDatabase.mat')
sce_kegg = kegg;
sce_swissprot = swissprot;
%Open Orthologs files
cd ../Annotation
fID        = fopen('AllGroupOrthologs.txt');
AllGroups  = textscan(fID,'%s %s %s %s %s','delimiter','\t');
fclose(fID);

fID        = fopen('singleCopyOrthologs.txt');
singleCopy = textscan(fID,'%s','delimiter','\n');
singleCopy = singleCopy{1};
n = length(singleCopy);
fclose(fID);
%Load Yali files for interconversion between W29 and CLIB122 IDs
fID        = fopen('Yali_annotation.txt');
yali = textscan(fID,'%s %s %s %s %s','delimiter','\t');
fclose(fID);
fID = fopen('yaliW29.txt');
w29 = textscan(fID,'%s %s','delimiter','\t');
fclose(fID);

%%Extract single-copy Ortho groups from AllGroups
OGids    = {};
genes    = cell(n,3);
uniprots = cell(n,3);
for i=1:n
    OGroup = singleCopy{i};
    index  = find(strcmpi(AllGroups{1},OGroup));
    if ~isempty(index)
        OGids      = [OGids; OGroup];
        matchStr   = AllGroups{5}{index};
        %Save Kmarx gene
        geneID     = AllGroups{3}(index);%findInProtDatabase(pID,kma_swissprot, kma_kegg);
        genes{i,2} = geneID;
        if ~isempty(matchStr)
            [genesCell, uniCell] = decomposeString(matchStr);
            uniprots{i,1} = uniCell{1};
            genes{i,1}    = genesCell{1};
            uniprots{i,3} = uniCell{3};
            genes{i,3}    = genesCell{3};
            %Save kma uniprots
            uniprots{i,2} = AllGroups{4}(index);
            %Find gene ID in database
        else
            matchStr  = AllGroups{4}{index};
            [genesCell,uniCell] = decomposeString(matchStr);
            uniprots{i,1} = uniCell{1};
            genes{i,1}    = genesCell{1};
            uniprots{i,3} = uniCell{3};
            genes{i,3}    = genesCell{3};
            %Save marxianus gene instead of its uniprot for missing cases
            uniprots{i,2} = AllGroups{3}(index);
       end
    end
    disp(['Ready with Orthogroup #' num2str(i)])
end
[genes,uniprots] = removeMissingValues(genes,uniprots,OGids);
cd ../Databases
genes = cell2table(genes);
writetable(genes,'orthologs_genes.txt','Delimiter','\t');
uniprots = cell2table(uniprots);
writetable(uniprots,'orthologs_proteins.txt','Delimiter','\t');
cd ../ComplementaryScripts
end

function [gene, pID] =  findInW29(pID,w29,yali)
index = find(strcmpi(w29{1},pID));
pID = '';
if ~isempty(index)
    index = index(1);
    gene = w29{2}{index};
    %Now translate the W29 id to CLIB122 gene ID
    if ~isempty(gene)
        index2 = find(strcmpi(yali{2},gene));
        if ~isempty(index2)
            %gene = yali{2}{index2};
            pID  = yali{4}{index2};
        end
    end
else
    gene = '';
end

end

function [genes, uniprots] = decomposeString(matchStr)
uniprots = cell(1,2);
genes    = cell(1,2);
matchStr = strsplit(matchStr,'|');
%Save sce uniprots
uniprots{1} = matchStr{2};
%Find gene ID in database
geneID     = findInProtDatabase(matchStr{2},sce_swissprot, sce_kegg);
genes{1} = geneID;

%Save yli uniprots
pID = matchStr{4};
%Find gene ID in database
[geneID,pID] = findInW29(pID,w29,yali);
genes{2}     = geneID;
uniprots{2}  = pID;
end

function gene = findInProtDatabase(pID,swissprot, kegg)
index = find(strcmpi(swissprot(:,1),pID));

if ~isempty(index)
    gene = swissprot{index,3};
    gene = strsplit(gene,' ');
    gene = gene{1};
else
    index = find(strcmpi(kegg(:,1),pID));
    if ~isempty(index)
        gene = kegg{index,3};
    else
        gene = '';
    end
end
end


function [filteredG,filteredP] = removeMissingValues(genes,proteins,OGids)
indexesProt = [];
indexesGen  = [];
for i=1:length(OGids)
    prot = proteins(i,:);
    gen  = genes(i,:);
    %Save those proteins with IDs for the three organisms
    if sum(~cellfun(@isempty,prot))==3
        indexesProt = [indexesProt; i];
    end
    %Save those genes with IDs for the three organisms
    if sum(~cellfun(@isempty,gen))==3
        indexesGen = [indexesGen; i];
    end
end
filteredG = [OGids(indexesGen), genes(indexesGen,:)];
filteredP = [OGids(indexesProt),proteins(indexesProt,:)];
end


