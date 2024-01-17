import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pymatgen.core.structure import Structure
from pymatgen.ext.matproj import MPRester
import numpy as np


####### Task A #######

df = pd.read_json("WF-CE_database_58332.json")
#print(df) #to check the extracted data

# Analyze distribution of column `cleavage_energy`  

sns.histplot(df['cleavage_energy'], kde = True)
plt.title('Distribution of Cleavage Energy')
plt.xlabel('Cleavage Energy')
plt.ylabel('Frequency')

#plot of cleavage energy as a pdf file
plt.savefig('cleavage_energy_distribution.pdf')
plt.show()

# Assess correlation between columns `WF` and `Fermi`
correlation = df[['WF', 'Fermi']].corr()

print("Correlation between WF and Fermi: ")
print(correlation)

# Ploting the correlation for visualization
sns.heatmap(correlation, annot=True, cmap="coolwarm")
plt.title("Correlation between WF and Fermi")
plt.savefig('Correlation_bt_WF_and_Fermi')
plt.show()


####### Task B #######

#Defining the function
def get_structure_info(struc: Structure):
    volume = struc.volume
    nelements = len(struc.composition.elements)
    
    # Calculate the angle between lattice vectors a and b 
    # (as there is no predefined parameter for that in pymatgen)
    a = struc.lattice.matrix[0]  # Lattice vector a
    b = struc.lattice.matrix[1]  # Lattice vector b
    cos_angle = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    angle_a_b = np.degrees(np.arccos(cos_angle))  # Convert to degrees

    return volume, nelements, angle_a_b

# Load mp-2490 crystal structure and feed it to the function `get_structure_info`

mpr = MPRester('wizRXXiTJ0TFSMkstmExr1y26b5DdOP0') #API key of mp-2490 structure
structure = mpr.get_structure_by_material_id("mp-2490") #making mp-2490 as pymatgen object 

# Testing the function with the loaded structure
volume, nelements, angle_a_b = get_structure_info(structure)
print(f"Volume: {volume} Å³")
print(f"Number of Unique Elements: {nelements}")
print(f"Angle between Lattice Vectors a and b: {angle_a_b} degrees")