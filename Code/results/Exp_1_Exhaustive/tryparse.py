from string import *
import os
import sys
import glob
exp_ids = [sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]] #Experiment id, Exp_1, Exp_2, Exp_3
fp_me=open(exp_ids[0]+"_ME.info","w")
fp_z = open(exp_ids[0]+"_Zrmse.info","w")
fp_fmax_infer = open(exp_ids[0]+"_Fmax_inferred.info","w")
fp_fmax_gd = open(exp_ids[0]+"_Fmax_gd.info","w")
fp_pr = open(exp_ids[0]+"_PR.info","w")
for i in range(1,len(exp_ids)):
    exp_id=exp_ids[i]
    path = "pspi_results_"+exp_id+"*.txt"
    files = glob.glob(path)
    files.sort()
    for input_file in files:
        fin = open(input_file,"r")
        count = 1
        for line in fin:
            if (count==52):
                outputstring = input_file+" "+line
                fp_me.write(outputstring)
            if (count==53):
                outputstring = input_file+" "+line
                fp_z.write(outputstring)
            if (count>53):
                if (count==54):
                    outputstring = input_file+" "+line
                    fp_fmax_infer.write(outputstring)
                if (count==55):
                    outputstring = input_file+" "+line
                    fp_fmax_gd.write(outputstring)
                if (count==58):
                    outputstring = input_file+" "+line
                    fp_pr.write(outputstring)
            count=count+1
fp_me.close()
fp_z.close()
fp_fmax_infer.close()
fp_fmax_gd.close()
fp_pr.close()

