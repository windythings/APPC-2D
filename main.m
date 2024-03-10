clear
close all

vecAng = 30;
alpha = 10;
CT = 1.5;

surfaceFiles{1} = sprintf('Airfoils/ONR-Coords/mainVec%d.dat',vecAng);
surfaceFiles{2} = sprintf('Airfoils/ONR-Coords/nacelleVec%d.dat',vecAng);
surfaceFiles{3} = 'Airfoils/ONR-Coords/krueger.dat';

for i = length(surfaceFiles):-1:1
    fileID = fopen(surfaceFiles{i});
    surfaces{i} = cell2mat(textscan(fileID,'%f%f','Delimiter',{'\t',','}));
    fclose(fileID);
end

options = wakeoptset(Wakesolver,'FunctionTolerance',1e-4,'Display','iter');
[Cl,Cd,Cp,xc] = Panel2D(surfaces,alpha,CT,options);

figure
hold on
for i = 1:length(Cp)
    plot(xc{i},Cp{i})
end
set(gca,'YDir','reverse')