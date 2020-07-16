function [ I,J ] = getcoordinate( index, Nx, Ny )

I=rem(index,Nx);

if rem(index,Nx)==0
    I=Nx;
end

J=(index-I)/Nx+1;

end

