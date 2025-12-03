# Charger la structure modifiée
load PAKP3_T0.pdb, genome

# Supprimer tout par défaut
hide everything, genome

# Représenter les segments comme sphères
show spheres, chain A
set sphere_scale, 0.2

#show spheres, chain B
#set sphere_scale, 0.4

# Représenter les liens CONECT comme sticks fins et transparents
show sticks, chain A
set stick_radius, 0.4
set stick_transparency, 0.3

#show sticks, chain B
#set stick_radius, 0.4
#set stick_transparency, 0.3

# coloriage du chromosome 
spectrum count, red_white_blue, chain A
#spectrum count, cyan_white_yellow, chain B

# make a selection with multiple indices
select specials, index 1+2+3+4+5+5613+5614+5615+5616+5617

# show them as larger black spheres
show spheres, specials
set sphere_scale, 1, specials
color red, specials

# label the origin
#select ori, index 1
#label ori, "ORI"
#set label_size, -10
#set label_placement_offset, [10.0, 10.0, -5.0]

#load other time points
load PAKP3_T3.pdb, genome
load PAKP3_T5.pdb, genome
load PAKP3_T7.pdb, genome
load PAKP3_T10.pdb, genome
load PAKP3_T13.pdb, genome
load PAKP3_T16.pdb, genome

# finition
bg_color white
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set ray_trace_disco, off
set depth_cue, off
set antialias, 2


