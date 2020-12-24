function z=compliance(x)
global number_of_nodes  nodes  nodal_force
global connectivity_matrix lengths_of_bars E fixed_nodes KG F U
%Initialisation à zero de la matrice de rigidité globale
KG = zeros(number_of_nodes*2);
%remplissage la matrice de rigidité globale 
pk=0;    
for i = 1: length(connectivity_matrix)
    for j = i+1: length(connectivity_matrix)
        if connectivity_matrix(i,j) == 1
            pk=pk+1;
            nx = double((nodes(j,1) - nodes(i,1))/lengths_of_bars(pk));
            ny = double((nodes(j,2) - nodes(i,2))/lengths_of_bars(pk));
            
            % matrice de rigidite locale
            K = [nx*nx,ny*nx,-nx*nx,-nx*ny;
                 ny*nx,ny*ny,-nx*ny,-ny*ny;
                 -nx*nx,-nx*ny,nx*nx,nx*ny;
                 -nx*ny,-ny*ny,nx*ny,ny*ny];
        
                 
            % section transversale de la pk ieme barre
            % lengths_of_bars(pk) longueur de la pk ieme barre
            
            K = ((E*x(pk))/lengths_of_bars(pk)^2)* K;
            
            
            
            %matrice de rigidité globale (assemblage)
            KG(2*i-1,2*i-1) = KG(2*i-1,2*i-1) + K(1,1);
            KG(2*i-1,2*i) = KG(2*i-1,2*i) + K(1,2);
            KG(2*i-1,2*j-1) = KG(2*i-1,2*j-1) + K(1,3);
            KG(2*i-1,2*j) = KG(2*i-1,2*j) + K(1,4);
            
            KG(2*i,2*i-1) = KG(2*i,2*i-1) + K(2,1);
            KG(2*i,2*i) = KG(2*i,2*i) + K(2,2);
            KG(2*i,2*j-1) = KG(2*i,2*j-1) + K(2,3);
            KG(2*i,2*j) = KG(2*i,2*j) + K(2,4);
            
            KG(2*j-1,2*i-1) = KG(2*j-1,2*i-1) + K(3,1);
            KG(2*j-1,2*i) = KG(2*j-1,2*i) + K(3,2);
            KG(2*j-1,2*j-1) = KG(2*j-1,2*j-1) + K(3,3);
            KG(2*j-1,2*j) = KG(2*j-1,2*j) + K(3,4);
            
            KG(2*j,2*i-1) = KG(2*j,2*i-1) + K(4,1);
            KG(2*j,2*i) = KG(2*j,2*i) + K(4,2);
            KG(2*j,2*j-1) = KG(2*j,2*j-1) + K(4,3);
            KG(2*j,2*j) = KG(2*j,2*j) + K(4,4);
        end
    end
end
% Vecteur chargement
F = zeros(2*number_of_nodes,1);
for i=1:number_of_nodes
    F(2*i-1)= nodal_force(i,1);
    F(2*i)= nodal_force(i,2);
end


for i = 1:length(KG)
    if (KG(i,:)==0) 
       KG(i,i)=1;
    end
end



%reduction de la dimension de la matrice de rigidité globale
KG(fixed_nodes,:)=[];
KG(:,fixed_nodes)=[];

%reduction de la dimension du vecteur chargement
F(fixed_nodes)=[];
warning('off','all')  % ! Pensez à traiter les configurations conduisant à une matrice singulière
C=rcond(KG);

if (C  > 10e-14)
U = KG\F;
%Energie de déformation
z = (transpose(F)*U)/2;
else
    z= 999999;
end;
