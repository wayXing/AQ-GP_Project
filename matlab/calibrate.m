function xCalibrated = calibrate(x,model)

if strcmp(model,'PMS5003')
  xCalibrated = 0.7778*x+2.6536;
%   xCalibrated = -67.0241*log(-0.00985*x+0.973658);
elseif strcmp(model,'PMS1003')
  xCalibrated = 0.5431*x+1.0607;
%   xCalibrated = -54.9149*log(-0.00765*x+0.981971);
elseif strcmp(model,'H1.1')
  xCalibrated = 0.4528*x+3.526;
end

end