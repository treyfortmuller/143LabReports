% Author: Andrew Jalbert
% This script is free to use
% Developed for EE 143 Lab Report 2 Fall 2018 at UC Berkeley

% Script solves for unkown Na given listed parameters below
% This will only work for a P type body of silicon
% Uses the MOSCAP model of a semiconductor device
syms x
% Charge of electron, measured in coulombs
q = 1.60217662E-19;
% permittivity of silicon in F/cm
e_si = 12*8.85418782E-14;
% Boltzmann constant times room temperature in Kelvin
% Should be units of joules
kT = 300 * 1.38064852E-23;
% intrinsic silicon, units of cm^-3
n_i = 1E10;

% Adjust the values below to meet your values
    % Voltage values in either eV or V
    Vt = 1;
    Vfb = -2;
    % oxide capacitance in F/cm^2
    Cox = 2.360E-8;

% Equation that will be solved for x
eqn = Vfb + 2*kT/q*log(x/n_i) + sqrt(2*q*2*e_si*x*kT/q*log(x/n_i))/Cox == Vt;

% Attempts to solve for x symbolically, but fails so it gives numerical
% value instead. The numerical value should be correct
solx = solve(eqn,x)
% If in doubt, plot the given function and find x such that equation is
% satisfied.


% One might have to divide solx by some power, like 1E16, to make it look
% more nice