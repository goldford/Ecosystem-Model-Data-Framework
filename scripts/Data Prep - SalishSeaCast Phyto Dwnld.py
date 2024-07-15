# Create Jul 11 2024
# Author: G Oldford
# Purpose: Download SalishSeaCast biology fields
# Notes:
#   - Based on ERDAPP site url generator https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.html
#   - NC files are approx 1 Gb / yr
import requests
from helpers import buildSortableString

# downloading the file using griddap
# note that download of 2007 (earliest year) didn't work
# yr_st = 2008
# yr_en = 2009
# mo_st = 1
# mo_en = 12
# date_st = str(yr_st) + "-" + buildSortableString(mo_st,1)
# date_en = str(yr_en) + "-" + buildSortableString(mo_en,1)
#
# # approx 1 Gb per year
# # URL generated from https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.html
# # https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.nc?ciliates[(2008-01-01T12:00:00Z):1:(2009-01-01T12:00:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)],diatoms[(2008-01-01T12:00:00Z):1:(2009-01-01T12:00:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)],flagellates[(2008-01-01T12:00:00Z):1:(2009-01-01T12:00:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)]
# griddap_url = ("https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.nc?" +
#                "ciliates[(" + date_st + "-01T12:00:00Z):1:(" + date_en + "-01T12:00:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)]," +
#                "diatoms[(" + date_st + "-01T12:00:00Z):1:(" + date_en + "-01T12:00:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)]," +
#                "flagellates[(" + date_st + "-01T12:00:00Z):1:(" + date_en + "-01T12:00:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)]"
#                )

# Define the range of years to download
start_year = 2007
end_year = 2019

for year in range(start_year, end_year):

    mo_st = 1
    da_st = 1
    mo_en = 12
    da_en = 31
    date_st = str(year) + "-" + buildSortableString(mo_st, 2) + "-" + buildSortableString(da_st,2)
    date_en = str(year) + "-" + buildSortableString(mo_en, 2) + "-" + buildSortableString(da_en,2)
    #date_st = str(yr_st) + "-" + str(mo_st)
    #date_en = str(yr_en) + "-" + str(mo_en)

    # Generate the URL
    griddap_url = (
        "https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.nc?"
        "ciliates[(" + date_st + "T12:00:00Z):1:(" + date_en + "T23:59:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)],"
        "diatoms[(" + date_st + "T12:00:00Z):1:(" + date_en + "T23:59:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)],"
        "flagellates[(" + date_st + "T12:00:00Z):1:(" + date_en + "T23:59:00Z)][(0.5000003):1:(20)][(0.0):1:(897.0)][(0.0):1:(397.0)]"
    )

    print(griddap_url)
    # Download the file
    response = requests.get(griddap_url)
    if response.status_code == 200:
        filename = f"SalishSeaCast_biology_{year}.nc"
        with open(filename, 'wb') as f:
            f.write(response.content)
        print(f"Downloaded: {filename}")
    else:
        print(f"Failed to download data for {year}. HTTP Status code: {response.status_code}")


