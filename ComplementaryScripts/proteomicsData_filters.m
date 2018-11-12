function [outputData,filtered] = proteomicsData_filters(rawData,DataIDs,genes,filename,samples,varFlag)
%
% proteomicsData_filters
%
% Function that Filters a proteomics dataset, the algorithm takes the  
% proteinsthat are present in at least two biological triplicates for all
% the experimental conditions. Zeros are assigned to all the missing values 
% in the original file. To run this script first, an EXCEL file with the 
% proteomics datasets must be loaded to the workspace.
%
%
%
% Last edited. Ivan Domenzain  2018-05-18
    current = pwd;
    [m,n]          = size(rawData);
    conditions     = n/samples;
    avgd_data      = cell(1,conditions);
    commonAvrgs    = {};%cell(1,conditions);
    avgd_ids       = cell(1,conditions);
    avgd_genes     = cell(1,conditions);
    filtered.data  = [];
    filtered.IDs   = {};
    filtered.genes = {};
    spreading      = cell(1,conditions);
    tolerance      = 2;
    % Main loop that parses each of the measured proteins
    for i=  1:m
        % Row position of the first experimental condition
        triplets_pos =(1:samples);
        % Loop that analyses each of the experimental conditions for the i-th
        % protein
        condPresence = 0;
        meanValsRow = zeros(1,conditions);
        for j=1:conditions
            datum = rawData(i,triplets_pos);
            nonZeros = find(datum>0);
%             if j==2
%                 tolerance = 1;
%             end
            % If the protein is present in at least two replicates then it
            % is considered as part of the filtered dataset
            
            %Also discard proteins with high variability
            %dataSpreading = (median(datum))/((max(datum)-min(datum)));
            %Calculate Relative standard deviation (RSD) for the measured
            %triplet of values
            dataSpreading = std(datum)/mean(datum);
            spreading{j}  = [spreading{j}; dataSpreading];
            if varFlag 
                Filter_condition = (sum(nonZeros)>=tolerance && ~isempty(DataIDs{i}) && dataSpreading<=1) && std(datum)>1E-6;
            else
                Filter_condition = (sum(nonZeros)>=2 && ~isempty(DataIDs{i}));
            end
            if Filter_condition
                condPresence   = condPresence + 1;
                value          = mean(datum);
                avgd_data{j}   = [avgd_data{j}; value];
                meanValsRow(j) = value;
                avgd_ids{j}    = [avgd_ids{j}; DataIDs(i)];
                avgd_genes{j} = [avgd_genes{j}; genes(i)]; 
            end
            triplets_pos = (1:samples)+((j)*samples);
        end
        % If a protein was measured (more than 1 time) for all the 
        % conditions, then the row is saved 
        if condPresence == conditions && ~isempty(DataIDs{i}) && ~isempty(genes{i}) %&& ~isequal(rawData(i,p(1)),rawData(i,p(2)),rawData(i,p(3)))
            filtered.data  = [filtered.data; rawData(i,:)];
            filtered.IDs   = [filtered.IDs; DataIDs(i)];
            filtered.genes = [filtered.genes; genes(i)];
            commonAvrgs    = [commonAvrgs; meanValsRow(:)];
        end
    end
    %Construct table with common measurements (all conditions) for
    %triplicates
    filtered.data = num2cell(filtered.data);
    %T             = horzcat(filtered.IDs,filtered.data);
    T             = horzcat(filtered.genes,filtered.data);
    
    T             = table(T);
    writetable(T,filename,'Delimiter','\t');
    %Construct table with common measurements (all conditions) for
    %mean value 
    T             = horzcat(filtered.genes,commonAvrgs);%
    T             = table(T);
    extensionPos  = strfind(filename,'.txt');
    
    filename = [filename(1:extensionPos-1) '_meanValues.txt'];
    writetable(T,filename,'Delimiter','\t');
    cd (current)
    %filtered.data = cell2mat(filtered.data);

    outputData.data  = avgd_data;
    outputData.ids   = avgd_ids;
    outputData.genes = avgd_genes;
    outputData.spreading = spreading;
end
