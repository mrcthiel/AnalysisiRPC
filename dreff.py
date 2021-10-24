import sys,os

targetDir = './outputs/'
sn = sys.argv[1]
wp = sys.argv[2]

#HV=['6800','7000','7100','7200','7400']
HV=[]
eff=[]

# check output directory
if not os.path.exists(targetDir):
    print("There is no outputs")
    os.exit()

# look for files with run number(sn)
for files in sorted(os.listdir(targetDir)):

    if sn in files:
        #print(files)
        rf = open(targetDir+files,'r')
        l = rf.readlines()

        for line in l:
            # look for the working point
            if "Efficiency" in line and wp in line and not "err" in line:
                print("Efficiency "+line.strip().split(": ",1)[1])
                eff.append(float(line.strip().split(": ",1)[1])*100.0)

            # How do we get an error?
            err = 0.01
            # OR
            #if "err" in line:
            #    err = line.strip().split(": ",1)[1]

# get list of HV as an input
while True:
    if len(HV) != len(eff):
        del HV[:]
        HVlist = raw_input("The number of HV elements is supposed to be "+str(len(eff))+". Please enter "+str(len(eff))+" elements of the list separated by space:\n>>")
        HV = HVlist.split()
        #print HV
    else: break

# write txt for eff_curve.C All blanks are "TAB ^I"
wf = open('eff_input_SN'+sn+'_'+wp+'.txt','w')
wf.write("HV 	Eff	Eff_err")
for i in range(len(eff)):
    wf.write("\n"+str(HV[i])+"	"+str(round(eff[i],1))+"	"+str(err))
wf.close()

# run the eff_curve with input
os.system("root -l -b -q eff_curve.C"+"'("+sn+","+"\""+wp+"\")'")
