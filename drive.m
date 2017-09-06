%engin

nGrids=8;

assumedSoln='constant';
figureID=11;
[error_phi0_n_constant, order_phi_nMinus1_constant]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_constant, order_phi_nMinus1_constant]=converger_FS(assumedSoln,nGrids,figureID)

assumedSoln='linear';
figureID=12;
[error_phi0_n_linear, order_phi_nMinus1_linear]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_linear, order_phi_nMinus1_linear]=converger_FS(assumedSoln,nGrids,figureID)

assumedSoln='quadratic';
figureID=13;
[error_phi0_n_quadratic, order_phi_nMinus1_quadratic]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_quadratic, order_phi_nMinus1_quadratic]=converger_FS(assumedSoln,nGrids,figureID)

assumedSoln='plus1Sqrt';
figureID=14;
[error_phi0_n_plus1Sqrt, order_phi_nMinus1_plus1Sqrt]=converger_FS(assumedSoln,nGrids,figureID)
[error_phi0_n_plus1Sqrt, order_phi_nMinus1_plus1Sqrt]=converger_FS(assumedSoln,nGrids,figureID)