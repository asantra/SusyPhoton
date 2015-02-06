import os, sys, re, glob


anafile  = open('SusyLimitMASTER.C')
ana_line = anafile.readlines()
anafile.close()

hfile    = open('SusyLimitMASTER.h')
hline    = hfile.readlines()
hfile.close()

mGlu  = 800
mLSP  = 25   



for mGlu in range(800,1601):
  if mGlu%50 == 0: 
    for mLSP in range(25,1576):
      if mLSP%25 == 0 and mLSP%50!=0 and mLSP < mGlu: 
        flag2 = False
        flag1 = False
        flag3 = False
        flag4 = False
        flag5 = False
        myfile  = open('SusyLimit.C','w')
        myfile2 = open('SusyLimit.h','w')
        
        for line2 in ana_line:
          if not flag1 and 'XXX' in line2:
	    strglu = str(mGlu)
	    strLSP = str(mLSP)
	    line2 = line2.replace('XXX', strglu)
	    line2 = line2.replace('YYY', strLSP)
	    #print "after replace: ",line2
	    flag1 = True
	  if not flag4 and 'GGG' in line2:
	    strglu = str(mGlu)
	    line2 = line2.replace('GGG', strglu)
	    flag4 = True 
	  if not flag5 and 'LLL' in line2:
	    strLSP = str(mLSP)
	    line2 = line2.replace('LLL', strLSP)
	    flag5 = True   
	  myfile.write(line2)
	myfile.close()
	
        for line3 in hline:
          if not flag2 and 'ZZZ' in line3:
            for line in open('SusySignal.txt'): 
	      susy = '/eos/uscms/store/user/asantra4/SusySignal/hist_analysis_SUSYSignal_'+str(mGlu)+'_'+str(mLSP)+'.root'
	      if susy in line:
		#print susy
		linea = line[:-1]
                line3 = line3.replace('ZZZ', linea)
            flag2 = True
          if not flag3 and 'TTT' in line3:
            for line in open('SusySignal.txt'): 
	      susy = '/eos/uscms/store/user/asantra4/SusySignal/hist_analysis_SUSYSignal_'+str(mGlu)+'_'+str(mLSP)+'.root'
	      if susy in line:
		lineb = line[:-1]
                line3 = line3.replace('TTT', lineb)
            flag3 = True
          myfile2.write(line3)
        myfile2.close()
            
           
          
        #myfile.close()
        os.system("bash SusyBash.sh")
        #os.system(".L SusyLimit.C++")
        #os.system("SusyLimit Ha")
        #os.system("Ha.Loop()")
        os.system("mv Signal_Gluino"+str(mGlu)+"_LSP"+str(mLSP)+"_AsymmetricPt.root /eos/uscms/store/user/asantra4/SusyComplete") # for asymmetricpt cut
        
 
    
      
     
      
	
	
    
      
    
 
   
#os.system('root -l -b -q '+'ana_'+str(ifile)+'.C')
#os.system('mv hist_analysis_SUSYSignal_'+str(ifile)+'.root /eos/uscms/store/user/asantra4/')
   
