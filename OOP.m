% OOP calculate orientational order parameter of a matrix
%
% Use this function to get orientational order parameter of actin stained
% images
%
% Inputs = nonzero_orientation (array of the nonzero orientation angles)
% Outputs = Orientation_order_parameter
function [ Orientation_order_parameter ] = OOP( nonzero_orientation )

%Tensor Method for Orientational Order Parameter
%Calculate x and y components of each vector r
if length(nonzero_orientation) > 2
    Vectors(1,:) = cos(nonzero_orientation);
    Vectors(2,:) = sin(nonzero_orientation);
    
    %Calculate the Orientational Order Tensor for each r and the average
    %Orientational Order Tensor (OOT_Mean)
    for i=1:2
        for j=1:2
            OOT_All(i,j,:)=Vectors(i,:).*Vectors(j,:);
            OOT_Mean(i,j) = mean(OOT_All(i,j,:));
        end
    end
    %Normalize the orientational Order Tensor (OOT), this is necessary to get the
    %order paramter in the range from 0 to 1
    OOT = 2.*OOT_Mean-eye(2);
    %Find the eigenvalues (orientational parameters) and
    %eigenvectors (directions) of the Orientational Order Tensor
    [~,orient_parameters]=eig(OOT);
    %orientational order parameters is the maximal eigenvalue, while the
    %direcotor is the corresponding eigenvector
    Orientation_order_parameter = max(max(orient_parameters));
else
    Orientation_order_parameter = 0;
end
end

