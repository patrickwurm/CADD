close all 
clear all
clc

additional_input ={};
additional_input.fe_damping_coefficient = 2.5e11;
additional_input.example_id_1 = 'Tensile_Test';
additional_input.example_id_2 = 'Test_1';
additional_input.fe_damping_coefficient = 10*2.5e11;
additional_input.selective_fe_damping = 1;

%Here, different random number seeds can be specified, to change the
%initial velocity distributions
%seeds = [1 263 9182 12020 13050 14555 15000 16777 17620 18845 19637];
seeds = 50000;

for i = 1:length(seeds)
disp(['Seed: ',num2str(seeds(i))])
additional_input.rng_seed = seeds(i);

additional_input.temperature = 50;
CADD2D('Tensile_Test_rotated','dynamic',additional_input)
CADD2D('Tensile_Test_rotated','static',additional_input)
CADD2D('Tensile_Test_rotated','hybrid',additional_input)

additional_input.temperature = 150;

CADD2D('Tensile_Test_rotated','dynamic',additional_input)
CADD2D('Tensile_Test_rotated','static',additional_input)
CADD2D('Tensile_Test_rotated','hybrid',additional_input)

additional_input.temperature = 250;

CADD2D('Tensile_Test_rotated','dynamic',additional_input)
CADD2D('Tensile_Test_rotated','static',additional_input)
CADD2D('Tensile_Test_rotated','hybrid',additional_input)
end
