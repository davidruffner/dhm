;batch file to record the details about this run
;This time we'll use hd5f files which seem like they will be
;very useful.

;Input information
path = './'
filename=path+'details.h5'

date = " 2014/07/24"

path2 = "../"
vidname = path2+'VTS_01_1.VOB';using the background as the video
bgvidname = path2+'VTS_01_1.VOB'
;dcvidname = path2+'VTS_22_1.VOB'
notes1 = ["Point trap stiffness ",$
        " 1.2um SiO2 spheres in DI water in deep cell",$
        " 447nm Cube holographic illumination",$
        " On Zappa with holoeye slm", $
        " 100x obj. 1x,  Henrique Moyses Data"]
ap = 0.6;particle size
np = 1.424;partice refractive index
lambda = .447;wavelength in vacuum
mpp=.135;microns per pixel
fps = 30;frames per secon

;temperature
Temp = 24.+273.;Read from computer with sensor by door to lab
kb = 1.3806488*10.d^(-23+18);pN*um/K ;Boltzmann constant
kT = kb*Temp ;pN*um; 

nm = refractiveindex(lambda,temp-273.) ;Medium refractive index



notes2 = ["using the background as the video VTS_01_1.VOB"]
;Start and end frames of video to analyze
startfrms = [0]
endfrms = [-1]
bg2use = [0];indgen(n_elements(startfrms)-1)
;Start and end frames of video for background
startfrmsbg = [0]
endfrmsbg = [-1]

;alpha 
alpha=1.0

;Range of interest
roi = [50,400,150,479]

;Set Power
power = 0.06;W
 

;Setup structure
struct = {date:date,$
          notes1:notes1,$
          notes2:notes2,$
          vidname:vidname,$
          bgvidname:bgvidname,$
          ;dcvidname:dcvidname,$
          ap:ap,$
          np:np,$
          lambda:lambda,$
          mpp:mpp,$
          fps:fps,$
          temp:temp,$
          kT:kT,$
          nm:nm,$
          startfrms:startfrms,$
          endfrms:endfrms,$
          bg2use:bg2use,$
          startfrmsbg:startfrmsbg,$
          endfrmsbg:endfrmsbg, $
          alpha:alpha, $
          roi:roi, $
          power:power $
          }

;Setup the hd5f file and write it
fileID = H5F_CREATE(filename)
datatypeID = H5t_idl_create(struct)
dataspaceid = h5s_create_simple(1)
datasetid = h5d_create(fileID,'mydata',datatypeid,dataspaceid)
h5d_write,datasetid,struct
h5f_close,fileid

 
