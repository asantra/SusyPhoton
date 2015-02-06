import os, sys, re, glob


anafile  = open('ana_MASTER.C')
ana_line = anafile.readlines()
anafile.close()


mGlu  = 1200
mLSP  = 25   


  
for mGlu in range(1200,1301):
  if mGlu%50 == 0: 
    for mLSP in range(25,1276):
      if mLSP%25 == 0 and mLSP%50!=0 and mLSP < mGlu: 
        flag2 = False
        flag1 = False
        myfile = open('ana_'+str(mGlu)+'_'+str(mLSP)+'.C','w')
        for line2 in ana_line:
          if not flag1 and 'ana_MASTER' in line2:
	    #print "before replace: ",line2
	    line2 = line2.replace('ana_MASTER', 'ana_'+str(mGlu)+'_'+str(mLSP))
	    #print "after replace: ",line2
	    flag1 = True
          if not flag2 and 'Check' in line2:
	    line2 = line2.replace('Check',str(mGlu)+'_'+str(mLSP))
            flag2 = True
          if '//X' in line2:
	    itext = 0
            for line in open('SUSYSignal-mGl1200to1300_mLSP-25to1275.txt'): 
	      susy = 'susyEvents_T5gg_'+str(mGlu)+'_'+str(mLSP)
	      if susy in line:
		#print susy
                infile = 'chain.Add(\"'+line
                infile2 = infile[:-1]+'\");'
                change = '//X'+str(itext)+'$'
                line2 = line2.replace(change, infile2)
                itext += 1
           
          myfile.write(line2)
        myfile.close()
        os.system("root -l -b -q ana_"+str(mGlu)+"_"+str(mLSP)+".C")
        os.system("mv hist_analysis_SUSYSignal_"+str(mGlu)+"_"+str(mLSP)+".root /eos/uscms/store/user/asantra4/SusySignal")
        os.system("rm ana_"+str(mGlu)+"_"+str(mLSP)+".C")
        
 
    
      
     
      
	
	
    
      
    
 
   
#os.system('root -l -b -q '+'ana_'+str(ifile)+'.C')
#os.system('mv hist_analysis_SUSYSignal_'+str(ifile)+'.root /eos/uscms/store/user/asantra4/')
   
