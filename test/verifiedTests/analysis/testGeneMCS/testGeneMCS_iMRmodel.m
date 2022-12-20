% The COBRAToolbox: testGeneMCS.m
%
% Purpose:
%     - Check the function of geneMCS using different solvers
%
% Authors:
%     - Luis V. Valcarcel 2022-12-12
%     - Naroa Barrena 2022-12-12

global CBTDIR
requiredSolvers = {'ibm_cplex'};
prepareTest('requiredSolvers',requiredSolvers);

% save the current path
currentDir = pwd;

% initialize the test
testDir = fileparts(which('testGeneMCS_iMRmodel'));
cd(testDir);

% define the solver packages to be used to run this test
solverPkgs = {'ibm_cplex'};

% Create Toy Example
ToyReactionIdentifiers = {'R1', 'R2', 'R3', 'R4', 'RBio'};
ToyReactionNames = {'Reaction 1', 'Reaction 2', 'Reaction 3', 'Reaction 4', 'Reaction Bio'};
ToyReactionFormulas = {' -> Met1', 'Met1 -> Met2', 'Met2 -> Met3', 'Met1 -> Met3', 'Met3 -> '};
ToyGrRuleList = {'g1', 'g2', 'g2', '(g2 and g3) or (g2 and g4) or g5', ''};

model = createModel(ToyReactionIdentifiers, ToyReactionNames, ToyReactionFormulas, 'grRuleList', ToyGrRuleList);
model = changeObjective(model, 'RBio');

% Create table with Boolean regulatory rules
ToySignalingNetwork = table({'g5', 'g7', 'g6', 'g7', 'g8', 'g7'}',...
    {'g2', 'g3', 'g4', 'g4', 'g5', 'g2'}', [+1 -1 +1 +1 +1 +1]',...
    'VariableNames', {'source_ENSEMBL', 'target_ENSEMBL', 'interaction'});


% expected solution
true_gmcs_0_layers = cell(2,1);
true_gmcs_0_layers{1} = {'g1'}';
true_gmcs_0_layers{2} = {'g2' 'g5'}';

true_gmcs_1_layers = cell(4,1);
true_gmcs_1_layers{1} = {'g1'}';
true_gmcs_1_layers{2} = {'g2' 'g5'}';
true_gmcs_1_layers{3} = {'g2' 'g8'}';
true_gmcs_1_layers{4} = {'g5' 'g7'}';


for k = 1:length(solverPkgs)
    fprintf(' -- Running testGeneMCS using the solver interface: %s ... ', solverPkgs{k});
    
    solverLPOK = changeCobraSolver(solverPkgs{k}, 'MILP', 0);
    
    if solverLPOK
        % Eliminate G-matrix if it exist
        if exist([testDir filesep 'G_toy_example_gMCS_iMRmodel_0_layers.mat'], 'file')
            delete G_toy_example_gMCS_iMRmodel_0_layers.mat
        end
        if exist([testDir filesep 'G_toy_example_gMCS_iMRmodel_1_layers.mat'], 'file')
            delete G_toy_example_gMCS_iMRmodel_1_layers.mat
        end
        if exist([testDir filesep 'CobraMILPSolver.log'], 'file')
            delete CobraMILPSolver.log
        end
        if exist([testDir filesep 'MILPProblem.mat'], 'file')
            delete MILPProblem.mat
        end
        if exist([testDir filesep 'tmp.mat'], 'file')
            delete tmp.mat
        end
        
        % Check errors when missing argument
        assert(verifyCobraFunctionError('calculateGeneMCS', 'inputs', {model, 20, 5}));
        assert(verifyCobraFunctionError('calculateGeneMCS', 'inputs', {'toy_example_gMCS', [], 20,5}));
        assert(verifyCobraFunctionError('calculateGeneMCS', 'inputs', {'toy_example_gMCS', model, [],5}));
        assert(verifyCobraFunctionError('calculateGeneMCS', 'inputs', {'toy_example_gMCS', model, 20}));
        
        % Calculate GMCS
        buildGmatrix_iMRmodel('toy_example_gMCS_iMRmodel', model, ToySignalingNetwork, 1);
        [gmcs, gmcs_time] = calculateGeneMCS('toy_example_gMCS_iMRmodel_1_layers', model, 20, 2);
        
        % Check if the solution obtained is the same as the expected
        % solution
        gmcsIsTrue = zeros(size(true_gmcs_1_layers));
        for i = 1:numel(gmcs)
            for j=1:numel(true_gmcs_1_layers)
                aux1 = gmcs{i};
                aux2 = true_gmcs_1_layers{j};
                if isequal(aux1,aux2)
                    gmcsIsTrue(j) = gmcsIsTrue(j)+1;
                    break
                end
            end
        end
        assert(sum(~logical(gmcsIsTrue))==0);
        % cellfun(@strjoin, gmcs, repmat({','}, numel(gmcs), 1), 'UniformOutput', 0)
        % cellfun(@strjoin, true_gmcs_0_layers, repmat({','}, numel(true_gmcs_0_layers), 1), 'UniformOutput', 0)
        % cellfun(@strjoin, true_gmcs_1_layers, repmat({','}, numel(true_gmcs_1_layers), 1), 'UniformOutput', 0)
        
        % Check the same for 0 layers (AKA, normal geneMCS)
        buildGmatrix_iMRmodel('toy_example_gMCS_iMRmodel', model, ToySignalingNetwork, 0);
        [gmcs, gmcs_time] = calculateGeneMCS('toy_example_gMCS_iMRmodel_0_layers', model, 20, 2);
        gmcsIsTrue = zeros(size(true_gmcs_0_layers));
        for i = 1:numel(gmcs)
            for j=1:numel(true_gmcs_0_layers)
                aux1 = gmcs{i};
                aux2 = true_gmcs_0_layers{j};
                if isequal(aux1,aux2)
                    gmcsIsTrue(j) = gmcsIsTrue(j)+1;
                    break
                end
            end
        end
        assert(sum(~logical(gmcsIsTrue))==0);
        
    else
        warning('The test testGeneMCS_iMRmodel cannot run using the solver interface: %s. The solver interface is not installed or not configured properly.\n', solverPkgs{k});
    end
    
    % Eliminate generated files
    if exist([testDir filesep 'G_toy_example_gMCS_iMRmodel_0_layers.mat'], 'file')
        delete G_toy_example_gMCS_iMRmodel_0_layers.mat
    end
    if exist([testDir filesep 'G_toy_example_gMCS_iMRmodel_1_layers.mat'], 'file')
        delete G_toy_example_gMCS_iMRmodel_1_layers.mat
    end
    if exist([testDir filesep 'CobraMILPSolver.log'], 'file')
        delete CobraMILPSolver.log
    end
    if exist([testDir filesep 'MILPProblem.mat'], 'file')
        delete MILPProblem.mat
    end
    if exist([testDir filesep 'tmp.mat'], 'file')
        delete tmp.mat
    end
    
    % output a success message
    fprintf('Done.\n');
end

% change the directory
cd(currentDir)
