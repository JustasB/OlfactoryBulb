function epsilon = kernel(t, t_j, tau_m_in, tau_s_in, tau_m, tau_s)
    % t time array
    % t_j spike time
    % tau_m_in  inhibition constant time
    % tau_s_in  inhibition constant time
    % tau_m excitation constant time
    % tau_s exitation constant time
    % epsilon output current form
    
    t_effective = (t - t_j);
    if (ne(tau_m_in, 0))
        epsilon_in=(exp(-t_effective/tau_s_in)-exp(-t_effective/tau_m_in)).*heaviside(t_effective);
        epsilon_in(find(isnan(epsilon_in)==1))=0;
        epsilon_in=(10/(1-tau_m_in/tau_s_in))*epsilon_in;
    else
        epsilon_in = zeros(1, length(t));
    end
    
    epsilon=(exp(-t_effective/tau_s)-exp(-t_effective/tau_m)).*heaviside(t_effective);
    epsilon(find(isnan(epsilon)==1))=0;
    epsilon=(10/(1-tau_m/tau_s))*epsilon;
    epsilon = epsilon-epsilon_in;
end