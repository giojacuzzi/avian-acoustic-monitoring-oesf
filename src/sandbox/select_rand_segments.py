## Randomly sample a number of files per species and copy to a new directory,
## given a root directory with audio segments from the classifier(s)
# 
# Input:
#
# - N number of samples to randomly select per species
# - Output directory 
# - Input root directory with directory tree assuming form below: 
#
# in_dir
#   Deployment
#       SiteUnit
#           SiteUnitDate
#               Species
#                   [.wav files]
#
in_dir = '/Volumes/gioj/OESF_processed/2022/segments'#'/Users/giojacuzzi/Desktop/OESF_processed_example/segments'
out_dir = '/Volumes/OESF_PAM_1/OESF_samples/2022' #'/Users/giojacuzzi/Downloads/elephant'
N = 500
##############################################################################
import os
import shutil
import pandas as pd

# Create output directory
os.makedirs(out_dir, exist_ok=True)

# Store file references per species
species_filepaths = []

# Iterate through directory tree
for deployment in os.listdir(in_dir):
    deployment_path = os.path.join(in_dir, deployment)
    if not os.path.isdir(deployment_path):
        continue
    print(f'Deployment {deployment}')
    
    for unit in sorted(os.listdir(deployment_path)):
        unit_path = os.path.join(deployment_path, unit)
        if not os.path.isdir(unit_path):
            continue
        print(f'    Unit {unit}')

        for date in sorted(os.listdir(unit_path)):
            date_path = os.path.join(unit_path, date)
            if not os.path.isdir(date_path):
                continue
            print(f'        Date {date}')

            for species in sorted(os.listdir(date_path)):
                species_path = os.path.join(date_path, species)
                if not os.path.isdir(species_path):
                    continue
                
                # Gather all audio files for the species
                files = [os.path.join(species_path, f) for f in os.listdir(species_path) if f.endswith('.wav')]
                print(f'            Species {species} ({len(files)})')

                # Store for later selection
                species_filepaths.extend([{'species': species, 'filepath': file} for file in files])

species_segments = pd.DataFrame(species_filepaths)
print("Finished collecting all segments:")
print(species_segments.head())

# Randomly sample files per species
print(f"Randomly sampling {N} segments per species...")
selected_files = species_segments.groupby('species', group_keys=False).apply(lambda x: x.sample(min(N, len(x)))).reset_index(drop=False)

# Create species index column
selected_files['species_idx'] = selected_files.groupby('species').cumcount() + 1

# Copy samples to output directory
print(f"Copying samples to output directory...")
for i, row in selected_files.iterrows():
    species_dest = os.path.join(out_dir, row['species'])
    os.makedirs(species_dest, exist_ok=True)
    out_file = f"{row['species_idx']}_{os.path.basename(row['filepath'])}"
    shutil.copy(row['filepath'], os.path.join(species_dest, out_file))

# Save selection reference table
selected_files['file'] = selected_files['filepath'].apply(os.path.basename)
selected_files[['species', 'file']].to_csv(os.path.join(out_dir, '_select_rand_segments.csv'), index=False)

print(f"Random selections copied to {out_dir}")