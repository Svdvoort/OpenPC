OpenPC is a toolbox for the creation and evaluation of Multi-Element Generalized Polynomial Chaos Expansions

QUICK START
----------------

To use OpenPC to make the PCE of a problem two things need to be provided: A settings file and a matlab function to interface with the 'blackbox'

Three example settings files are provided in the settings folder.
Example_settings.json contains the settings for a simple 4 dimensional problem. 
As reference for the meaning of the different terms it is suggested to read: Perkó, Zoltán, et al. "Fast and accurate sensitivity analysis of IMPT treatment plans using Polynomial Chaos Expansion." Physics in medicine and biology 61.12 (2016): 4646.
Example_settings_arbitrary.json contain the settings for a problem where an arbitrary PDF is defined, which can be found in Multi/PDF/
Example_settings_ME.json for the construction of a multi-element Polynomial Chaos Expansion.

Each of the settings has an example run in Examples/, first run initOpenPC from the main directory and then run one of the examples.

The corresponding blackbox interfaces for the settings and the example runs can be found in Blackbox/Functions

NOTATION
----------------

The code has some notation specifications:
- functions, variables start with a lower case letter, structures start with a upper case letter.
- There are a few shorthands used for easier reading:
	- N: Number of
	- u: unique
	- I: index of
	- LI: logical index of
	- i: i-th element of (used in iterations)
  For example the variable N_dims means: Number of dimensions. 
  These shorthands can also be combined, e.g. N_u_pol_types means the number of unique polynomial types.




