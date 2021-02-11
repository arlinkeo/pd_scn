cd 'C:/Users/Arlin/surfdrive/pd_imaging_scn/GeneCategoryEnrichmentAnalysis-0.1.1'

% Defining null phenotype ensemble
enrichmentParams = GiveMeDefaultEnrichmentParams();
enrichmentParams.fileNameOut
enrichmentParams.dataSource = 'human-direct';

% AHBA entrez IDs
geneInfo = readtable('../../AHBA_Arlin/probe_info_2018-11-18.csv');
ahba_entrezID = geneInfo{:, 'entrez_id'};% get entrez IDs
ahba_entrezID = num2cell(ahba_entrezID);
ahba_entrezID = cell2mat(ahba_entrezID);
ahba_entrezID = ahba_entrezID';

% load AHBA expr data
donors = ["donor9861", "donor10021", "donor12876", "donor14380", "donor15496", "donor15697"];
m = [];
for i = 1:length(donors)
    d = donors{i};
    disp(d);
    f = strcat("../../AHBA_Arlin/gene_expr_", d, "_2018-11-18.csv");
    e = readtable(f, 'ReadVariableNames', true, 'ReadRowNames',true);%, 'VariableNamingRule', 'preserve');
    e.Properties.VariableNames = strcat(d, "_", e.Properties.VariableNames);
    m = [m e];
end
m = table2cell(m);
m = cell2mat(m);
m = m';
s = struct('expressionMatrix', m); % region x genes
s.entrezIDs = ahba_entrezID;

% Load phenotype data (binary)
p = readtable('../pd_scn/output/network_phenotype.csv', 'ReadVariableNames', true, 'ReadRowNames',true);
networkC = p{ "Network_C", :}';
networkD = p{ "Network_D", :}';

% Compute null distribution
ComputeAllCategoryNulls(s,enrichmentParams,[],true,true);

% Compute enrichment of phenotype
GOTablePhenotype = EnsembleEnrichment(m, '../GCEA_results/Network_C',networkC);
GOTablePhenotype = EnsembleEnrichment(m, '../GCEA_results/Network_D',networkD);







GOTable = SingleEnrichment(geneScores,ahba_entrezID,enrichmentSettings);