%% Instruction
  % to run different cases, change the manufacturer only!
%% Info
% Converger is for grid refinement analysis. 
% Right now, it needs to be aware of the geometry and the material. So it
% can call the manufacturer to get the MMS problem and reference solution. 
% The geometry and material also need to be passed to the physics modules, 
% because it's part of the problem description. 
% It needs to know the geometry and is responsible for generating the grid
% and pass the grid information to the coupler. 
% It also calculates the error and reveal the rate of convergence. 
function [order_phi0]=...
  converger_MoC_LS_Eig(assumedSoln,k_MMS,nGrids,refinementRatio,N,angErrorRemoval)
format long;
% Case configure options
if ~exist('assumedSoln','var')
  assumedSoln='const-const';
  assumedSoln='flat-expMu';
  assumedSoln='linear-expMu';
  assumedSoln='quadratic-expMu';
  assumedSoln='plus1Sqrt-expMu';
%   assumedSoln='sine-complex';
%   assumedSoln='nonSeparableDivision';
%   assumedSoln='nonSeparableAddition';
end
if ~exist('k_MMS','var')
  k_MMS=1.02;
end
if ~exist('nGrids','var')
  nGrids=8%8%4%4%6;%10;%8;
end
if ~exist('refinementRatio','var')
  refinementRatio=2;
end
if ~exist('N','var')
  N=4; % angular discretization 
end
if ~exist('angErrorRemoval','var')
  angErrorRemoval='no';
  angErrorRemoval='partial';
  angErrorRemoval='complete';
end

% Geometry
Tau=10; 

gridMeshSize_iGrid=zeros(nGrids,1);
error_phi0_iGrid=zeros(nGrids,1);
error_k_iGrid=zeros(nGrids,1);

for iGrid=1:nGrids
  J=5*refinementRatio^iGrid;
%   J=20;
  gridMeshSize_iGrid(iGrid)=Tau/J;
  iGrid
  % Material
  field1='Sig_t_j';          value1=ones(J,1);
  field2='Sig_ss_j';         value2=ones(J,1)*0.4;
  field3='Sig_gamma_j';      value3=ones(J,1)*0.5;
  field4='Sig_f_j';          value4=ones(J,1)*0.1;
  field5='nuSig_f_j';        value5=ones(J,1)*0.2;
  field6='thermal_cond_k_j'; value6=ones(J,1);
  field7='kappaSig_f_j';     value7=ones(J,1)*0.1; % kappa=1.0;
  mat = struct(field1,value1,field2,value2,field3,value3,... 
    field4,value4,field5,value5,field6,value6,field7,value7);

  [phi0_MMS_j,...
    psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n,...
    error_ang_j,error_hat_ang_j,...
    phi0_guess_j,k_guess]=... 
        manufacturer_MoC_LS_Eig(J,N,Tau,mat,assumedSoln,k_MMS);
  
  % Call eigen solver
  if strcmp(angErrorRemoval,'no')
    error_ang_j=error_ang_j*0.0;
    error_hat_ang_j=error_hat_ang_j*0.0;
    [phi0_j,k]=MoC_LS_Eig_module(J,N,Tau,mat,...
      psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n,...
      error_ang_j,error_hat_ang_j,...
      phi0_guess_j,k_guess);
  end
  if strcmp(angErrorRemoval,'partial')
    [phi0_j,k]=MoC_LS_Eig_module(J,N,Tau,mat,...
      psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n,...
      error_ang_j*0.0,error_hat_ang_j*0.0,...
      phi0_guess_j,k_guess);
  end
  if strcmp(angErrorRemoval,'complete')
    [phi0_j,k]=MoC_LS_Eig_module(J,N,Tau,mat,...
      psi_b1_n,psi_b2_n,Q_MMS_j_n,Q_MMS_hat_j_n,...
      error_ang_j,error_hat_ang_j,...
      phi0_guess_j,k_guess);
  end

  error_phi0_iGrid(iGrid)=norm(phi0_j-phi0_MMS_j,2)/sqrt(J)
  error_k_iGrid(iGrid)=k...
    -k_MMS
  
end

% Calculate the order of accuracy
order_phi0_nMinus1=ones(nGrids-1,1);
order_k_nMinus1=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi0_nMinus1(j)=log(error_phi0_iGrid(j)/error_phi0_iGrid(j+1)) / ...
    log(gridMeshSize_iGrid(j)/gridMeshSize_iGrid(j+1));
  order_k_nMinus1(j)=log(error_k_iGrid(j)/error_k_iGrid(j+1)) / ...
    log(gridMeshSize_iGrid(j)/gridMeshSize_iGrid(j+1));
end

%% Visualize the asymptotic convergence
orderPlotGrid=[gridMeshSize_iGrid(1) gridMeshSize_iGrid(end)];

scalarFluxErrorRMS_plot_handle=figure;
loglog(gridMeshSize_iGrid,error_phi0_iGrid,'*');
title({'cell-averaged scalar flux error convergence',...
  ['\phi_{MMS}: ' assumedSoln '; k_{MMS}: ' num2str(k_MMS)]});
xlabel('mesh size [cm]');
ylabel('RMS error of scalar flux');

hold on;
orderGuess=round(order_phi0_nMinus1(end));
errorStt=error_phi0_iGrid(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'r--');
loglog(orderPlotGrid,secondOrder,'g--');
loglog(orderPlotGrid,thirdOrder,'b--');
loglog(orderPlotGrid,fourthOrder,'k--');
legend('scalar flux error','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');

set(get(gca,'xlabel'),'FontName','Times New Roman');
set(get(gca,'ylabel'),'FontName','Times New Roman');
set(get(gca,'title'),'FontName','Times New Roman');
set(findobj(gcf, 'Type', 'Legend'),'FontName','Times New Roman');

hold off;
savefig(scalarFluxErrorRMS_plot_handle,['MoC_LS_Eig_Physor2018_phi_complete_AER_' assumedSoln]);
% savefig(scalarFluxErrorRMS_plot_handle,['MoCEig_Physor2016_cellAveraged_phi_convergence_' assumedSoln]);

kError_plot_handle=figure;
loglog(gridMeshSize_iGrid,error_k_iGrid,'*');
title({'k error convergence',...
  ['\phi_{MMS}: ' assumedSoln '; k_{MMS}: ' num2str(k_MMS)]});
xlabel('mesh size [cm]');
ylabel('k error');

hold on;
orderGuess=round(order_k_nMinus1(end));
errorStt=error_k_iGrid(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'r--');
loglog(orderPlotGrid,secondOrder,'g--');
loglog(orderPlotGrid,thirdOrder,'b--');
loglog(orderPlotGrid,fourthOrder,'k--');
legend('k error','1st Order','2nd Order',...;
  '3rd Order','4th Order','location','best');

set(get(gca,'xlabel'),'FontName','Times New Roman');
set(get(gca,'ylabel'),'FontName','Times New Roman');
set(get(gca,'title'),'FontName','Times New Roman');
set(findobj(gcf, 'Type', 'Legend'),'FontName','Times New Roman');

hold off;
savefig(kError_plot_handle,['MoC_LS_Eig_Physor2018_k_complete_AER_' assumedSoln]);
% savefig(kError_plot_handle,['MoCEig_Physor2016_k_convergence_' assumedSoln]);

%% Dispaly the result
% Display the problem description and results
disp '=================';
display(['assumedSoln: ' assumedSoln]);
display(['Number of grids: ' num2str(nGrids)]);
display(['refinementRatio: ' num2str(refinementRatio)]);
display(['quad set order: ' num2str(N)]);

error_phi0_iGrid
order_phi0_nMinus1
error_k_iGrid
order_k_nMinus1

order_phi0=order_phi0_nMinus1(end);
order_k=order_k_nMinus1(end);
