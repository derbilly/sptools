function newData = fInterpolate(data, frequency, newFrequency)
magnitude = abs(data);
phase = unwrap(angle(data));
newMagnitude = interp1(frequency,magnitude,newFrequency,'pchip'); 
newPhase = interp1(frequency,phase,newFrequency,'pchip');
newData = newMagnitude.*exp(1j*newPhase);
end
