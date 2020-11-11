function [V,D] =Hamiltonian_QCA(matrix)
sigmaZ = [1,0;0,-1];
sigmaX = [0,1;1,0];
gamma =0.01;

% Define the massive kink energy matrix
cell_ind=find(matrix==5);
[m,n]=size(matrix);
num_cells=length(cell_ind);

[row,col] = ind2sub([m,n],cell_ind);
Kink_Matrix=zeros(length(cell_ind));

for i = 1:num_cells-1
        for j = i+1:num_cells
            xy1=[row(i),col(i)];
            xy2=[row(j), col(j)];
            Kink_Matrix(i,j)=get_Ek_scale((xy1(1)-xy2(1)),(xy1(2)-xy2(2)));
        end
end

    %%
% Create the first part of Hamiltonian, consisting of the tunneling
% barriers
H_gamma = 0;
for i = 1:num_cells
    H_gamma_new = -gamma*tensor_product(i,num_cells,sigmaX);
    H_gamma = H_gamma + H_gamma_new;
end

%%
% Create second part of Hamiltonian consisting of driver-cell interaction
driver_ind=find((matrix==1)|(matrix==-1));
[rowD,colD] = ind2sub([m,n],driver_ind);
xy1D= [rowD(1) colD(1)];
xy2D= [rowD(2) colD(2)];
xy3D= [rowD(3) colD(3)];

    H_Driver = 0;
    
    for i = 1:num_cells
        xy2=[row(i) col(i)];
        P1Ek= get_Ek_scale((xy1D(1)-xy2(1)),(xy1D(2)-xy2(2)))*matrix(rowD(1),colD(1))   ;
        P2Ek= get_Ek_scale((xy2D(1)-xy2(1)),(xy2D(2)-xy2(2)))*matrix(rowD(2),colD(2))   ;
        P3Ek= get_Ek_scale((xy3D(1)-xy2(1)),(xy3D(2)-xy2(2)))*matrix(rowD(3),colD(3))   ;
        driverEk = P1Ek + P2Ek + P3Ek;
        H_Driver_new = -0.5*driverEk*tensor_product(i,num_cells,sigmaZ);
        H_Driver = H_Driver + H_Driver_new;

    end
%% Create the third part of Hamiltonian, consisting of cell-cell coupling

H_coupl = 0;
    for i = 1:num_cells-1
        for j = i+1:num_cells
            H_coupl_new = -0.5*Kink_Matrix(i,j)*tensor_product(i,num_cells,sigmaZ)*tensor_product(j,num_cells,sigmaZ);
            H_coupl = H_coupl + H_coupl_new;
        end
    end
    

H = H_gamma+ H_Driver + H_coupl;

[V,D] = eig(H);
D = diag(D);
    
    
end 