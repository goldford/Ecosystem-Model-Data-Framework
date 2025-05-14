import os
import shutil

# paths
#source_folder = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//forcing//ECOSPACE_in_3day_vars_1980-2018//vartemp1_C_0-10mAvg"  #
# source_folder = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//forcing//ECOSPACE_in_3day_vars_1980-2018//varsalt1_PSU_0-10mAvg"  #
source_folder = "C://Users//Greig//Documents//GitHub//Ecosystem-Model-Data-Framework//data//forcing//NEMO_prepped_as_ASC_monthly//TempVertMean4m"
#source_folder = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//LTL_MODELS//RDRS forcings//Wind_RDRS//Ecospace//stress"

destination_folder = source_folder
os.makedirs(destination_folder, exist_ok=True)
files = os.listdir(source_folder)
# variable = "PAR-VarZ-VarK"
#variable = "vartemp1_C_0-10mAvg"
variable = "varsalt1_PSU_0-10mAvg"
#variable = "RDRS_windstress10m"


# Loop through each file and copy with the new name format
for file_name in files:


    if file_name.startswith(variable + "_1980") and file_name.endswith(".asc"):
        print(file_name)
        # Extract the day of the year (e.g., 02, 05, etc.)
        if variable == "vartemp1_C_0-10mAvg":
            day_of_year = file_name.split("_")[4].split(".")[0]
            print(day_of_year)
        elif variable == "varsalt1_PSU_0-10mAvg":
            day_of_year = file_name.split("_")[4].split(".")[0]
            print(day_of_year)
        elif variable == "RDRS_windstress10m":
            day_of_year = file_name.split("_")[3].split(".")[0]
            print(day_of_year)
        else:
            day_of_year = file_name.split("_")[2].split(".")[0]

        # Construct the new file name for 1979
        new_file_name = f"{variable}_1978_{day_of_year}_butreally_1980_{day_of_year}.asc"

        # Define full file paths for source and destination
        source_file_path = os.path.join(source_folder, file_name)
        destination_file_path = os.path.join(destination_folder, new_file_name)

        # Copy the file to the destination folder with the new name
        shutil.copy(source_file_path, destination_file_path)

print("Script executed successfully!")