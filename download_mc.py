import urllib.request

url = "https://raw.githubusercontent.com/smistad/marching-cubes/master/mc_tables.h"
try:
    response = urllib.request.urlopen(url)
    data = response.read().decode('utf-8')
    with open("SurfaceDeterminationProject/mc_tables.h", "w") as f:
        f.write(data)
    print("mc_tables.h downloaded successfully.")
except Exception as e:
    print(f"Failed to download: {e}")
