function bool = free_energy_improved(free_energy, new_free_energy, options)
    diff = new_free_energy - free_energy;
    if abs(diff/free_energy) < options.threshold %change is small
        bool = 0;
    elseif diff < 0  %free energy is worse
        bool = 0;
    elseif diff == 0 %free energy is same
        bool = 0;
    else  % free energy is significantly improved
        bool = 1;
    end
end