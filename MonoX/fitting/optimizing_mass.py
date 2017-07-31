import sys
import os



#mass=[("50","75","0"),("75","100","1"),("100","150","2"),("150","3000","3")]
mass=[("100","150","2")]
PAT=os.getenv('PWD')
#os.system("source "+PAT+"/setup.sh")
for i,j,k in mass:
	os.system("cp "+PAT+"/../python/Selection_doubleb_XX.py "+PAT+"/../python/Selection_doubleb.py")
	os.system("sed -i 's/XX_MS1_XX/"+i+"/g' "+PAT+"/../python/Selection_doubleb.py")
	os.system("sed -i 's/XX_MS2_XX/"+j+"/g' "+PAT+"/../python/Selection_doubleb.py")
	os.system("mkdir -p /data/t3home000/mcremone/lpc/jorgem/skim/monohiggs_boosted/limits_mass"+k+"/")
	os.system("rm /data/t3home000/mcremone/lpc/jorgem/skim/monohiggs_boosted/limits_mass"+k+"/*.root")
	os.system("cp makeLimitForest_monohiggs_boosted_XX.py makeLimitForest_monohiggs_boosted.py")
	os.system("sed -i 's/XX_FD_XX/limits_mass"+k+"/g' makeLimitForest_monohiggs_boosted_n2ddt.py")
	os.system("cp makeLimitForest_XX.sh makeLimitForest.sh")
	os.system("sed -i 's/XX_LIM_XX/limits_mass"+k+"/g' makeLimitForest.sh")
	os.system("sed -i 's/XX_MS_XX/"+k+"/g' makeLimitForest.sh")
	os.system("source "+PAT+"/makeLimitForest.sh boosted")
	os.system("echo limits_mass"+k+" is done")
