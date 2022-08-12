% -- Description --
% This is a simple script to load the current simulated using the TDDFT code
% by Simon Madsen.  It calculates the nonlinear current generated within a 
% chain of atoms.  Here we are processing the data from all atomic orbitals
% to determine the net generated current and determin the high-harmonic
% radiation from it.

%-- Settings -- 
N_o = 200; %no. of orbitals
N_t = 95000; %no. of time steps 
dt = 0.1;
savefile = 'totalcurrent.dat';

%Read the whole matrix in
t = [0:N_t-1]*dt;

orbital_currents = dlmread('currentsup.dat', ' ');
orbital_currents = reshape(orbital_currents(:, 2), N_o, N_t);

totalcurrent = sum(orbital_currents);

data_out = [t', totalcurrent'];
csvwrite(savefile, data_out);

figure();
plot(t, totalcurrent);




