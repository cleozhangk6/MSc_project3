
# pymol input.pdb -d orient -c input.pml -g input.png
#   -c : read commands from input file
#   -d : orient molecule for 'best' view
#   -g : output image file

# Display selection
hide everything, all
show cartoon, all

# Background colour
bg_color white

# Colour our residues
spectrum b, red_blue

# Image Size
viewport 800, 800

#  Output image settings
set depth_cue,1
set antialias, 1
set cartoon_fancy_sheets, 1
set two_sided_lighting, on
set cartoon_loop_quality=100
set ray_opaque_background, on
set ray_shadows,off

# Raytrace image for best quality
ray


