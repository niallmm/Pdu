% extract data out of folders
% originally programmed to extract a jc, kc sweep and plot contour plot of
% both p and a

% add paths
add_paths
changeplot

% the base folder where your results are saved 
basefolder = '/Users/niallmangan/Dropbox/CCMtesting/';
% folder name where plots will be stored
foldername = ['000plotfolder'];
% make a folder to save plots
mkdir([basefolder foldername])

% find all the folders in this path
A = dir(basefolder);

for i = 1:length(A)
    
    foldername = A(i).name
    

    
    % if results are in the folder load results, load parameters and calculate things 
        % deal with Mac directory structure nonsense
    if foldername(1) == '.'
        foldername
    elseif exist([basefolder foldername '/res.mat'])>0
        foldername
         % load results all referenced by res.resultname
         C = load([basefolder foldername '/res.mat']);
         res = C.res;
         % load parameters
         B = load([basefolder foldername '/p.mat']);
         p = B.p;
         
         % create vectors of things you are interested in plotting
         % we will re-organize these vectors later
         
         jc_vec(i) = p.jc;
         kc_vec(i) = p.kcA;
         
         p_MCP_mM_vec(i) = res.p_MCP_mM;
         a_MCP_mM_vec(i) = res.a_MCP_mM;
         
     end
end
         


% take the vectors, find unique values of your varied parameters, then find
% values of resulst at these varied parameters and put in matrix
% points do not need to be regularly spaced
[jc, kc, pM] = xyzcontour(jc_vec',kc_vec',p_MCP_mM_vec');   
    
a = figure;
contourf(jc, kc, pM,20, 'ShowText', 'off', 'LineWidth', 2)
% can also use 'contour' to plot without fillled in color
xlabel('active transport rate of diol [cm/s]')
ylabel('permeability of compartment [cm/s]')
title('Concentration of p')
set(gca,'yscale','log');
set(gca,'xscale','log');
colorbar
colormapnew % change color scheme
savefig([basefolder foldername '/jcvskc_p.fig'])
saveas(a, [basefolder foldername '/jcvskc_p.png'], 'png')  


% take the vectors, find unique values of your varied parameters, then find
% values of resulst at these varied parameters and put in matrix
[jc, kc, aM] = xyzcontour(jc_vec',kc_vec',a_MCP_mM_vec');   
    
b = figure;
contourf(jc, kc, aM,20, 'ShowText', 'off', 'LineWidth', 2)
xlabel('active transport rate of diol [cm/s]')
ylabel('permeability of compartment [cm/s]')
title('Concentration of p')
set(gca,'yscale','log');
set(gca,'xscale','log');
colorbar
colormapnew % change color scheme
savefig([basefolder foldername '/jcvskc_a.fig'])
saveas(a, [basefolder foldername '/jcvskc_a.png'], 'png')    