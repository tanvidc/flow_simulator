function [ index ] = getindex( I,J,Nx,Ny )
roworder=1; columnorder=0;

if roworder==1
    index=(J-1)*Nx+I;
end

if columnorder==1
    index= (I-1)*Ny+J;
end

end

