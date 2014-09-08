;Get background for first part of video

;;;;;Input parameters from Spheretool;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Region of interest
rc = [ 371, 341]
;rad =          80
;;; Flags
fixap =           0 
fixzp =           0
fixnp =           0
fixnm =           1
fixalpha =           1
;;; Starting parameters for fits
;zp =       150.
path = './'
in = h5_parse(path+'details.h5')
details = in.mydata._data

ap =      details.ap
np =      details.np 
alpha =   details.alpha     
lambda =  details.lambda
temp =    details.temp
nm =      details.nm
mpp =     details.mpp
roi =     details.roi

vidname =   details.vidname
bgvidname = details.bgvidname
dcvidname = details.dcvidname
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;Combining Data into one folder of frames;;;;;;;
w=656
h=480
;Check if the video is ripped
vob = (strsplit(vidname,"_",/extract))[1]
f = file_search("../vid"+vob+"_*frames"+"/*/*/*.gdf")
;if not then rip video
if f[0] eq '' then begin $
   print,"ripping video.."  & $
   vidnames = file_search(vidname)  & $
   vob = (strsplit(vidname,"_",/extract))[1] & $
   fname="../vid"+vob+"_1frames" & $
   foreach vidname,vidnames do begin $
      print,"ripping "+vidname & $
      parts = strsplit(vidname,"_.",/extract) & $
      fname="../vid"+parts[2]+"_"+parts[3]+"frames" & $
      mp_ripper,vidname,fname,dimensions = [w,h], /greyscale & $
   endforeach & $
endif



;;;;;;;;;;;;;;;;;;;;get clean image of first frame
f = file_search("../vid"+vob+"_*frames"+"/*/*/*.gdf")
a = (read_gdf(f[1]))[5:644,*]
write_image,"../frame1.png","png",a

;;;;;;;;;;;;;;;;;;;;;;;Getting backgrounds
;Check if the bg video is ripped
bgvob = (strsplit(bgvidname,"_",/extract))[1]
bgf = file_search("../vid"+bgvob+"_*bgframes"+"/*/*/*.gdf")
;if not then rip bg video
if bgf[0] eq '' then begin $
   bgfname="../vid"+bgvob+"_1bgframes" & $
   mp_ripper,bgvidname,bgfname,dimensions = [w,h], /greyscale & $
   bgf = file_search("../vid"+bgvob+"_*bgframes"+"/*/*/*.gdf") & $
endif
;bgf = f  ;Un-comment to use the video as the background

startfrmsbg = details.startfrmsbg
endfrmsbg =   details.endfrmsbg
nsteps = n_elements(startfrmsbg)
print,"getting backgrounds"
for i=0,nsteps-1 do begin $
   if (file_search("../data/background"+vob+"part*.gdf"))[0] $
      ne '' then begin $
         print,"already got background" & $
         break & $
   endif & $
   print,"start",startfrmsbg[i] & $
   print,"end",endfrmsbg[i] & $
   ;Background frames
   
   getbackground,frames=bgf,$
              fout="../data/background"+vob+"part"+strtrim(i,2)+".gdf",$
              startfrm=startfrmsbg[i],endframe=endfrmsbg[i],step=5 & $
   frame = (read_gdf(f[startfrmsbg[i]]))[5:644,*] & $
   bg = read_gdf("../data/background"+vob+"part"+strtrim(i,2)+".gdf") & $
   help,frame  & $ 
   help,bg & $
   frameclean = (frame-32.)/float(bg-32. >0.0000000001) & $
   plotimage,bytscl(frameclean),/iso & $
   write_gdf, frameclean,$
             "../data/frame"+strtrim(endfrmsbg[i],2)+$
                       "clean"+vob+"part"+strtrim(i,2)+".gdf" & $
   write_image,"../figures/frame"+strtrim(endfrmsbg[i],2)+$
                       "clean"+vob+"part"+strtrim(i,2)+".png","png",$
   bytscl(frameclean) & $
   write_gdf, frame,$
             "../data/frame"+strtrim(endfrmsbg[i],2)+$
                       "raw"+vob+"part"+strtrim(i,2)+".gdf" & $
   write_image,"../figures/frame"+strtrim(endfrmsbg[i],2)+$
                       "raw"+vob+"part"+strtrim(i,2)+".png","png",$
   bytscl(frame) & $
   print,"finished getting background "+ strtrim(i,2) & $
   endfor
print,"moving on"

;;;;;;;;;;;;;;;;;;;;;;;;Get dark count field
w=656
h=480
;Check if the dc video is ripped
dcvob = (strsplit(dcvidname,"_",/extract))[1]
dcf = file_search("../vid"+dcvob+"_*dcframes"+"/*/*/*.gdf")
help,dcf
if dcf[0] eq '' then begin $
   print,"ripping video for dark counts.."  & $
   vidnames = file_search(dcvidname)  & $
   foreach vidname,vidnames do begin $
      print,"ripping "+vidname & $
      parts = strsplit(vidname,"_.",/extract) & $
      fname="../vid"+parts[2]+"_"+parts[3]+"dcframes" & $
      mp_ripper,vidname,fname,dimensions = [w,h], /greyscale & $
   endforeach & $
   dcf = file_search("../vid"+dcvob+"_*dcframes"+"/*/*/*.gdf") & $
endif
dc = (read_gdf(dcf[1]))[5:644,*] 
write_image,"../dcframe.png","png",dc 

;; ;;;;;;;;;;;;;;;;;;;;;;Analyzing data
featurethreshold=100
zp = 108.
rad=zp*tan(35*!pi/180.)<180
pickn=1
tag = ""
startfrms = details.startfrms
endfrms =   details.endfrms
bg2use =    details.bg2use
nparts = n_elements(startfrms)
for i=0,nparts-1 do begin $ 
   bgfile = "../data/background"+vob+"part"+strtrim(bg2use[i],2)+".gdf" & $
   fout = "../data/streamdata"+vob+"part"+strtrim(i,2) & $
   startfrm=startfrms[i] & $
   endframe=endfrms[i] & $
   p1=dhmstream_mp(frames=f,$
                bgfile=bgfile,$
                fout=fout, $
                lambda=lambda,mpp=mpp,$
                nm = nm,np=np,ap=ap,zp=zp,alpha=alpha,rad=rad,$
                roi=roi,featurethreshold=featurethreshold,$
                startfrm = startfrm, endframe = endframe,$
                fixap=0,fixzp=0,fixnp=0,fixalpha=1,$
                pickn=pickn,dc=dc) & $
endfor





