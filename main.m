% This code calculates the polarization of a cell of choice in a QCA array.
% In this case, we have implemented a majority gate. The physical layout is
% described in 'matrix', with 1 or -1 representing driver cells, 0
% representing empty cells, and x as a placeholder for normal cells. 

% The matrix cells are indexed as [row, col] in the order of increasing
% column index values. For the majority gate circuit simulated in this
% code, this means: 

%         0 0 0 0 0 D  0
%         0 0 0 0 0 5  0
%         0 0 0 0 0 6  0
%         0 0 0 0 0 7  0
%         0 0 0 0 0 8  0
%         D 1 2 3 4 9 11
%         0 0 0 0 0 10 0
%         0 0 0 0 0 D  0
% If we want to find the resulting polarization of the output cell, we must
% input index 11 into the polarization calculator. 
%
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start of Code
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=5;              % arbitrary placeholder for the cell location in the matrix
output_state=11;  % index for the output state whose polarization we want to know

matrix=[0 0 0 0 0 -1 0; ...
        0 0 0 0 0 x 0; ...
        0 0 0 0 0 x 0; ...
        0 0 0 0 0 x 0; ...
        0 0 0 0 0 x 0; ...
        1 x x x x x x; ...
        0 0 0 0 0 x 0; ...
        0 0 0 0 0 -1 0];
    
cell_ind=find(matrix==x);   % number of cells in the simulation 
[m,n]=size(matrix) ;
sigmaZ = [1,0;0,-1];

[V,D] = Hamiltonian_QCA(matrix); % call the function that calculates the total Hamiltonian 

% find the eigenvector that corresponds to the minimum eigenvalue (ground
% state)
min_index= find(min(D));
ground_eigenV= V(:,min_index);

% calculate polarization of desired output state by finding the expectation
% value of the PauliZ matrix for that state.
polarization= (ground_eigenV')*tensor_product(output_state,length(cell_ind),sigmaZ)*ground_eigenV;
disp(polarization)
