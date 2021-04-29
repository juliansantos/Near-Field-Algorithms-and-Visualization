% here enter the reduced matrix from zero towards chess

function prop_matrix = calculatePropagationMatrix(x_mesh, y_mesh, layers, lambda)

lambda_mm = lambda*1e3;
k = x_mesh.^2+y_mesh.^2;
delta_z = sqrt(k+layers(2)^2) - sqrt(k+layers(1)^2);
prop_matrix = exp(-1i*2*pi*delta_z/lambda_mm);

end 