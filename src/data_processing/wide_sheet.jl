using CSV
using DataFrames

# Define the base directory
cd(@__DIR__)
base_dir = "../../result/occ"

# Define the 8 variants (folder names)
variants = [
    "employment", "hourly_rate", "hours", "income", 
    "income_share", "inequality", "median", "unemployment"
]

# Iterate through each variant folder
for var in variants
    folder_path = joinpath(base_dir, var)
    
    # Check if the directory exists
    if !isdir(folder_path)
        println("Warning: Directory not found: $folder_path. Skipping...")
        continue
    end
    
    # Get all CSV files in the directory
    files = filter(f -> endswith(f, ".csv") && startswith(f, "irf_"), readdir(folder_path))
    
    if isempty(files)
        println("Warning: No CSV files found in $folder_path. Skipping...")
        continue
    end
    
    # Read the first file to initialize the base DataFrame with the horizon column
    first_file_path = joinpath(folder_path, files[1])
    first_df = CSV.read(first_file_path, DataFrame)
    
    # Create an initial wide DataFrame containing only the horizon column
    wide_df = DataFrame(horizon = first_df.horizon)
    
    # Iterate through all occupation CSV files for this variant
    for file in files
        file_path = joinpath(folder_path, file)
        
        # Extract occupation name from the filename 
        # (e.g., "irf_group1_Managerial.csv" becomes "group1_Managerial")
        occ_name = replace(file, r"^irf_" => "", r"\.csv$" => "")
        
        # Read the current occupation's CSV
        df = CSV.read(file_path, DataFrame)
        
        # Keep only horizon and beta_abs, renaming beta_abs to the occupation name
        temp_df = select(df, :horizon, :beta_abs => Symbol(occ_name))
        
        # Left join to append the data by horizon
        wide_df = leftjoin(wide_df, temp_df, on = :horizon)
    end
    
    # Sort by horizon to ensure correct chronological order
    sort!(wide_df, :horizon)
    
    # Generate the output file path
    out_file = joinpath(base_dir, "wide_$(var).csv")
    
    # Save the wide table as a CSV
    CSV.write(out_file, wide_df)
    println("Success: Generated wide table at $out_file")
end

println("Done: All 8 wide tables have been successfully generated!")