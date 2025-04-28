#!/bin/bash
#SBATCH -J temp_zarr2nc           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 20:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o logs/analysis.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e logs/analysis.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

conda activate zaar2netcdf

# Input and output paths
input_path="/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/Cape_Hatteras/spatial_long/dr_200"
output_path="/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/Cape_Hatteras/netcdf_files/dr_200"

# Ensure the output directory exists
mkdir -p "$output_path"

# Loop through all .zarr files in the input directory
for zarr_file in "$input_path"/*.zarr; do
    # Extract the base name of the file (without the directory)
    file_name=$(basename "$zarr_file")
    
    # Generate the output .nc file name
    nc_file="$output_path/${file_name%.zarr}.nc"
    
    # Convert the .zarr file to .nc
    echo "Converting $zarr_file to $nc_file"
    zarr2netcdf "$zarr_file" --output "$nc_file" --use_dask --verbose
done

# Input and output paths
input_path="/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/Cape_Hatteras/temporal_long/W_20"
output_path="/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/Cape_Hatteras/netcdf_files/W_20"

# Ensure the output directory exists
mkdir -p "$output_path"

# Loop through all .zarr files in the input directory
for zarr_file in "$input_path"/*.zarr; do
    # Extract the base name of the file (without the directory)
    file_name=$(basename "$zarr_file")
    
    # Generate the output .nc file name
    nc_file="$output_path/${file_name%.zarr}.nc"
    
    # Convert the .zarr file to .nc
    echo "Converting $zarr_file to $nc_file"
    zarr2netcdf "$zarr_file" --output "$nc_file" --use_dask --verbose
done
