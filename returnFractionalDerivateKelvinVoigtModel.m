function [modelFit, modelInput] = returnFractionalDerivateKelvinVoigtModel(Hz, modelFrequencyRange,currentStorageMod,currentLossMod)

model = @(x,xdata)(x(1)^(1-x(2)))*(1i*2*pi*xdata).^(x(2)) + x(3) + 1i.*x(4);

[mu,alpha,mu_z, mu_imag,resnorm,residual] = frac_derivative_KV_model(Hz,currentStorageMod,currentLossMod);
modelInput(1) = mu;
modelInput(2) = alpha;
modelInput(3) = mu_z;
modelInput(4) = mu_imag;
modelFit = model(modelInput, modelFrequencyRange);