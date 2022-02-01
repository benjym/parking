# This script will generate all of the figures required for the manuscript ``Parallel parking vehicle alignment strategies''

# Spatiotemporal figures (one for each simulation and then a collated plot)
python parking.py json/spatiotemporal.json
python make_comparison_spatiotemporal_plot.py json/spatiotemporal.json

# Strategies figure
python parking.py json/strategies.json
python summarise_strategies.py

# Road length figure
python parking.py json/blockface_lengths.json
python summarise_blockface_lengths.py

# SUVs figure (not included in paper)
python parking.py json/SUVs.json
python summarise_SUVs.py

# motorcycles figure (not included in paper)
python parking.py json/motorcycles.json
python summarise_motorcycles.py

# convert images to JPEG as per Findings requirements
# mogrify -density 400 -quality 90 -format jpg strategies.pdf blockface_lengths.pdf spatiotemporal_all.pdf
