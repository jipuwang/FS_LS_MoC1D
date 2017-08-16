%% Instruction
  % to run different cases, change the manufacturer only!
%% Info
% Grid refiner is for grid refinement analysis. 
% Right now, it needs to be aware of the geometry and the material. So it
% can call the manufacturer to get the MMS problem and solution
% The geometry and material also need to be passed to the coupler, so the
% coupler can keep passing the info on to the modules, because it's part of
% the problem description. 
% It needs to know the geometry and is responsible for generating the grid
% and pass the grid information to the coupler. 

clear;
nGrids=4%6;%10;%8;
refinementRatio=2;

% Geometry
Tau=10; 

% Case configure options
fbType='noFeedback'; % options: 'noFeedback','linear','squareRootPlus1'
assumedSoln='sine_sine'; % options: 'const_quadratic','sine_sine'

error_phi0_n=zeros(nGrids,1);
error_T_n=zeros(nGrids,1);
gridMeshSize_n=ones(nGrids,1);
N=8; % angular discretization, fixed not refined. 
for iGrid=1:nGrids
  J=5*refinementRatio^iGrid;
  gridMeshSize_n(iGrid)=Tau/J;
  
  
  %calculate the analytical solution for sqrt problem
  h=Tau/J;
  phi0_ana_j=ones(J,1);
  phi0_MMS =@(x) 2*sqrt(x+1);
  for j=1:J
    x_L=(j-1)*h;x_R=j*h;
    phi0_ana_j(j)=1/h*integral(phi0_MMS,x_L,x_R);
  end
  % call the coupler to solve the above manufactured problem
  [phi0_j]=MoC1D_prob5_LS_sqrt(J,N);
  
  % Calculate the error compared to manufactured solution
  error_phi0_n(iGrid)=norm(phi0_j-phi0_ana_j,2)/sqrt(J);

  %% Plot the solution
  x=linspace(0,10,J);
  scalarFlux_plot_handle=figure(13);
  clf;
  plot(x,phi0_j,'-');
  hold on;
  fplot(@(x) 2*sqrt(x+1), [0 10]);
  % title('scalar flux');
  xlabel('mesh size [cm]');
  ylabel('scalar flux');
  legend('numerical','analytical');
  hold off;

end

% Calculate the order of accuracy
order_phi_nMinus1=ones(nGrids-1,1);
for j=1:nGrids-1
  order_phi_nMinus1(j)=log(error_phi0_n(j)/error_phi0_n(j+1)) / ...
    log(gridMeshSize_n(j)/gridMeshSize_n(j+1));
end

% Display the result
error_phi0_n
order_phi_nMinus1


% Visualize the results
orderPlotGrid=[gridMeshSize_n(1) gridMeshSize_n(end)];

scalarFluxErrorRMS_plot_handle=figure(11);
clf;
loglog(gridMeshSize_n,error_phi0_n,'*');
% title('scalar flux error convergence');
xlabel('mesh size [cm]');
ylabel('scalar flux error RMS');

hold on;
orderGuess=round(order_phi_nMinus1(end));
errorStt=error_phi0_n(end)*refinementRatio^(orderGuess*(nGrids-1));
firstOrder=[errorStt errorStt/refinementRatio^(nGrids-1)];
secondOrder=[errorStt errorStt/refinementRatio^(2*(nGrids-1))];
thirdOrder=[errorStt errorStt/refinementRatio^(3*(nGrids-1))];
fourthOrder=[errorStt errorStt/refinementRatio^(4*(nGrids-1))];
loglog(orderPlotGrid,firstOrder,'--');
loglog(orderPlotGrid,secondOrder,'--');
loglog(orderPlotGrid,thirdOrder,'--');
loglog(orderPlotGrid,fourthOrder,'--');
legend('scalar flux error','1st Order','2nd Order',...
  '3rd Order','4th Order','location','best');
hold off;

