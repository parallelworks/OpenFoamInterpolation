import numpy as np
from scipy.interpolate import griddata

class OpenFoamFields():
     def __init__(self):
       # Starting and ending times in the controlDict file
       self.startTime=0
       self.endTime=0
       # Time interval in time units between saved fields
       self.writeInt=0
       # Path to the OpenFOAM case
       self.path2case='def'
       # List with all the saved fields
       self.savedList=[]
       # Number of dimensions of the problem: 2 or 3.
       self.nD=3
     # Uses the path to the OpenFoam case to initialize
     # the class atributes. The user must also specify 
     # the number of dimensions of the problem: nD. 
     def initialize(self,path,nD):
       self.path2case=path
       self.nD=nD
       path2controlDict=self.path2case+"/system/controlDict"
       fcontrolDict=open(path2controlDict,'r')
       for line in fcontrolDict:
          words=line.split(' ');
          if words[0]=="startTime":
             startTime=''
             for c in line:
                if c.isdigit() or c=='.':
                  startTime=startTime+c
             startTime=float(startTime)
          elif words[0]=="endTime":
             endTime=''
             for c in line:
                if c.isdigit() or c=='.':
                  endTime=endTime+c
             endTime=float(endTime)
          elif words[0]=="deltaT":
             deltaT=''
             for c in line:
                if c.isdigit() or c=='.':
                  deltaT=deltaT+c
             deltaT=float(deltaT)
          elif words[0]=="writeInterval":
             writeInterval=''
             for c in line:
                if c.isdigit() or c=='.':
                  writeInterval=writeInterval+c
             writeInterval=int(writeInterval)
       self.writeInt=writeInterval*deltaT
       self.startTime=startTime
       self.endTime=endTime
       Nfields=int((self.endTime-self.startTime)/self.writeInt)
       self.savedList=np.zeros(Nfields+1,dtype=np.double)
       for k in range(Nfields+1): 
           self.savedList[k]=self.writeInt*k
  
     # Reads the instantaneous pressure and velocity fields saved under path2instant
     # into the array Fields. Note that OpenFoam stores the value of these fields in
     # the center of the cells. Therefore, the x, y and z coordinates of a any cell
     # e are stored in Fields[e,0], Fields[e,1] and Fields[e,2], respectively.
     # provide the x, y and z coordinates of the center of cell "e", respectively. 
     # Similarly, the Ux, Uy, Uz and P fields corresponding to cell "e" are stored in
     # Fields[e,3], Fields[e,4], Fields[e,5] and Fields[e,6]. 
     def getInstFields(self,path2instant):
       fccx=open(path2instant+'/ccx','r')
       Lfccx=fccx.readlines()
       fccy=open(path2instant+'/ccy','r')
       Lfccy=fccy.readlines()
       fccz=open(path2instant+'/ccz','r')
       Lfccz=fccz.readlines()
       fU=open(path2instant+'/U','r')
       LfU=fU.readlines()
       fp=open(path2instant+'/p','r')
       Lfp=fp.readlines()
       k=0
       startLine=0
       initField=0
       elem=0
       for line in LfU:
          if initField==0:
             words = line.split(' ');
             for word in words:
                if word=="internalField":
                  initField=1
                  startLine=k
          elif initField==1:
            nel=int(LfU[k])
            Fields=np.zeros((nel,7),dtype=np.double)
            initField=2
          elif initField==2 and k>startLine+2 and k<=startLine+2+nel:
            ccx=float(Lfccx[k])
            Fields[elem,0]=ccx
            ccy=float(Lfccy[k])
            Fields[elem,1]=ccy
            if self.nD==2:
               ccz=0.0
            elif self.nD==3:
               ccz=float(Lfccz[k])
            else:
               print("Wrong number of dimensions nD. Must be 2 or 3!!!")
            Fields[elem,2]=ccz
            line=line[line.find("(")+1:line.find(")")]
            velComponents=line.split(' ')
            n=3
            for velComponent in velComponents:
             Fields[elem,n]=float(velComponent)
             n=n+1
            p=float(Lfp[k])
            Fields[elem,6]=p
            elem=elem+1 
          k=k+1
       return Fields

     # Combines timeInterp with spatialInterp to provide the velocity and pressure fields at
     # any time t and point (x,y,z).
     def interpolate(self,t,x,y,z):
       # Interpolation in time:
       interpFields=self.timeInterp(t)
       # Interpolation in space
       txyz2UVWP=self.spatialInterp(x,y,x,interpFields)
       return txyz2UVWP

     # Interpolates between the saved fields to calculate their value at a given instant t.
     # If t>timespan of simulation, last saved field is extrapolated (steady state)
     # Bug:  Won't work with t between init conditions and first saved fields
     def timeInterp(self,t):
       # Interpolation in time:
       div_t=t/self.writeInt
       mod_t=int(div_t)
       if mod_t>0:
         if mod_t==div_t and mod_t<=self.endTime: # No interpolation is needed
           fieldName=str(self.savedList[mod_t])
           nameParts=fieldName.split('.')
           if len(nameParts)>1 and int(nameParts[1])==0:
             fieldName=nameParts[0]
           path2fields=self.path2case+"/"+fieldName
           interpFields=self.getInstFields(path2fields)
         # If time is beyond simulation timespan the last fields will be assumed
         # as steady state.
         elif mod_t>=self.endTime: # Extrapolation to latest field
           fieldName=str(max(self.savedList))
           nameParts=fieldName.split('.')
           if len(nameParts)>1 and int(nameParts[1])==0:
             fieldName=nameParts[0]
           path2fields=self.path2case+"/"+fieldName
           interpFields=self.getInstFields(path2fields)
         else: # Interpolation is needed
           alpha=(t/self.writeInt)-mod_t
           # Previous fields
           fieldName=str(self.savedList[mod_t])
           nameParts=fieldName.split('.')
           if len(nameParts)>1 and int(nameParts[1])==0:
             fieldName=nameParts[0]
           path2prev=self.path2case+"/"+fieldName
           prevFields=self.getInstFields(path2prev)
           # Next fields
           fieldName=str(self.savedList[mod_t+1])
           nameParts=fieldName.split('.')
           if len(nameParts)>1 and int(nameParts[1])==0:
             fieldName=nameParts[0]
           path2next=self.path2case+"/"+fieldName
           nextFields=self.getInstFields(path2next)
           interpFields=(1-alpha)*prevFields+alpha*nextFields
       else:
         print('Cannot use initial conditions to interpolate. Select a later time.')
       return interpFields
     # Calculates the value of the pressure and velicity fields at a given point (x,y,z)
     # by interpolating the value of these fields at the center of the cells (Fields).
     # If (x,y,z) is outside of the domain, the value nan is returned.
     def spatialInterp(self,x,y,z,Fields):
       # Interpolation in space
       txyz2UVWP=np.zeros(4,dtype=np.double)
       if self.nD==2:
         txyz2UVWP[0]=griddata(Fields[:,[0,1]],Fields[:,3],(x,y),method='linear')
         txyz2UVWP[1]=griddata(Fields[:,[0,1]],Fields[:,4],(x,y),method='linear')
         txyz2UVWP[2]=griddata(Fields[:,[0,1]],Fields[:,5],(x,y),method='linear')
         txyz2UVWP[3]=griddata(Fields[:,[0,1]],Fields[:,6],(x,y),method='linear')
       elif self.nD==3:
         txyz2UVWP=np.zeros(4,dtype=np.double)
         txyz2UVWP[0]=griddata(Fields[:,[0,1,2]],Fields[:,3],(x,y,z),method='linear')
         txyz2UVWP[1]=griddata(Fields[:,[0,1,2]],Fields[:,4],(x,y,z),method='linear')
         txyz2UVWP[2]=griddata(Fields[:,[0,1,2]],Fields[:,5],(x,y,z),method='linear')
         txyz2UVWP[3]=griddata(Fields[:,[0,1,2]],Fields[:,6],(x,y,z),method='linear')
       return txyz2UVWP
     # Provides the minimum MinAndMax[0] and maximum MinAndMax[1] of the desired field at time "t".
     # Available fields: P, Ux, Uy, Uz, Umag. 
     def  minmax(self,field,t):
       MinAndMax=np.zeros(2,dtype=np.double)
       # Interpolation in time:
       interpFields=self.timeInterp(t)
       if field=='P':
         MinAndMax[0]=np.amin(interpFields[:,6])
         MinAndMax[1]=np.amax(interpFields[:,6])
       elif field=='Ux':
         MinAndMax[0]=np.amin(interpFields[:,3])
         MinAndMax[1]=np.amax(interpFields[:,3])
       elif field=='Uy':
         MinAndMax[0]=np.amin(interpFields[:,4])
         MinAndMax[1]=np.amax(interpFields[:,4])
       elif field=='Uz':
         MinAndMax[0]=np.amin(interpFields[:,5])
         MinAndMax[1]=np.amax(interpFields[:,5])
       elif field=='Umag':
         nel=np.size(interpFields,0)
         Umag=np.zeros(nel,dtype=np.double)
         for e in range(nel):
            Umag[e]=np.sqrt(np.power(interpFields[e,3],2)+np.power(interpFields[e,4],2)+np.power(interpFields[e,5],2))
         MinAndMax[0]=np.amin(Umag)
         MinAndMax[1]=np.amax(Umag)
       else:
         print(("There is no field named: %s. Choose: Umag, Ux, Uy, Uz or P") % field)
       return MinAndMax
     # Provides the minimum MinAndMax[0] and maximum MinAndMax[1] of the desired field throught the timespan of the simulation.
     # Available fields: P, Ux, Uy, Uz and Umag.
     def allTimeMinmax(self,field):
        k=0
        for savedFields in self.savedList:
           if k==1:
             MinAndMax=self.minmax(field,savedFields)
           elif k>1:
             MinAndMax_new=self.minmax(field,savedFields)
             if MinAndMax_new[0]<MinAndMax[0]:
               MinAndMax[0]=MinAndMax_new[0]
             if MinAndMax_new[1]>MinAndMax[1]:
               MinAndMax[1]=MinAndMax_new[1]
           k=k+1
        return MinAndMax
