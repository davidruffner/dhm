;@dhmstream
;Purpose: Batch file for analyzing a video file with DHM. 

;Description: Takes a VOB file, along with input parameters, then
;calculates a background if it is not already created. Next it checks
;for a crash. The main part of the program is a loop through all the
;frames where in each frame features are found and fit with DHM.

;General Usage
;;run once to get a background
;; Use the background to clean an image and run spheretool to get
;input parameters, save dhmstream3 with those parameters in the folder
;with the videofile.
;Then run
 ;;;      IDL> @dhmstream

;;Output
    ; writes to file streamdatap1full.gdf
    ; p = read_gdf("streamdatap1full.gdf")

    ;p[0,0:2,*] - x,y,z positions
    ;p[0,3,*]   - radius, ap
    ;p[0,4,*]   - particle refractive index, np
    ;p[0,5,*]   - alpha
    ;p[0,6,*]   - medium refractive index, nm
    
    ;p[1,0:6,*] - errors
    
    ;p[0:1,7,*] - frame number, odd frames on the zero side i think
                           ; and even on the 1 side.

;;Input (inputs are typed directly into the program)
      ;  filename for video, ex. "VTS_01_1.VOB"
      ; starting values for the variables.

;;Requirements: 
      ; 
 
;;Version History:
    ;4_21_11 : Cleaned and documented DBR
    ;6_27_12 : stop using idlvideo DBR
    ;7_10_14 : DBR. Many things have changed over past years
               ;today I added the dc keyword for inputing dark count image 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function dhmstream_mp,path=path,bgfile=bgfile,fout=fout,$
                       lambda=lambda,mpp=mpp,$
                       nm = nm,np=np,ap=ap,zp=zp,alpha=alpha,rad=rad,$
                       roi=roi,featurethreshold=featurethreshold,$
                       startfrm = startfrm, endframe = endframe,$
                       pickn=pickn,frames=frames,$
                       fixalpha=fixalpha,fixdelta=fixdelta,$
                       fixap=fixap,fixzp=fixzp,fixnp=fixnp,$
                       randompos=randompos,calcfraction=calcfraction,$
                       dc=dc
;;;;;;;;;Set Initial Parameters;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;.r gpuinit
;gpuinit
gpu =0                  ;Use gpu acceleration when set
doshow =1       
if n_elements(path) eq 0 then path = "./vid01frames"   ;INPUT VIDEOFile  
if n_elements(bgfile) eq 0 then bgfile = "background.gdf"
if n_elements(fout) eq 0 then fout = "streamdata"  
f = file_search(path+"/*/*/*.gdf")
if n_elements(frames) ne 0 then f = frames 
nframes = n_elements(f)
;;;;;Input parameters from Spheretool;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Region of interest
rc = [ 371, 341]
if n_elements(rad) eq 0 then rad = 40
;;; Flags
if n_elements(fixap) eq 0 then fixap =           0 
;if n_elements(fixzp) eq 0 then fixzp =           0
if n_elements(fixnp) eq 0 then fixnp =           0
fixnm =           1
if n_elements(fixalpha) eq 0 then fixalpha =           1
if n_elements(fixdelta) eq 0 then fixdelta =           1
;;; Starting parameters for fits
if n_elements(zp) eq 0 then zp =       62.
if n_elements(ap) eq 0 then ap =      .75
if n_elements(np) eq 0 then np =      1.45
if n_elements(alpha) eq 0 then alpha =       1.0000000
if n_elements(lambda) eq 0 then lambda =      0.447
if n_elements(nm) eq 0 then $
       nm =       refractiveindex(lambda,24.);roomtemperature was 24C on 5-31
if n_elements(mpp) eq 0 then mpp =      0.13500000 
if n_elements(dc) eq 0 then dc = 32.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;nm = refractiveindex(lambda, 26.) ; refractive index of water
print,"nm=",nm

; region of interest
if n_elements(roi) eq 0 then roi = [0,639,0,479]
x0 = roi[0]
x1 = roi[1]
y0 = roi[2]
y1 = roi[3]
;; x0 = 0
;; x1 = 330
;; y0 = 50
;; y1 = 400

nx = x1-x0+1
ny = y1-y0+1

nmed =300                       ; number of frames in running median
buf = bytarr(nmed, nx, ny)  & $  ; frame buffer

if n_elements(featurethreshold) eq 0 then featurethreshold =100

if n_elements(startfrm) eq 0 then startfrm = 1
if n_elements(endframe) eq 0 then endframe = nframes-1
if endframe eq -1 then endframe = nframes-1

;;;;;;;;;;;;Checking for Crash;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print,"checking for crash"

g = file_search(fout+"?.gdf", count = oldfiles)
if oldfiles eq 2 then begin $ ;;set so it doesn't recheck
   print,"recovering" & $
   p = read_gdf(g[0]) & $
   t = read_gdf(g[1]) & $
   zp = p[0,2,-1] & $ ;Careful!! it estimates the radius from zp
   rad = zp*(100/244.)+10 & $ ; only works for one particle in field of view
   startfrm = total(t[0:1,-1]) +1  & $
   print,"startframe=",startfrm & $
endif else begin $
   p = fltarr(2, 9) & $
   t = fltarr(2) & $
   ;n0 = nmed & $
   print,"no crash" & $
   endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;Get Background;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
g = file_search(bgfile, count = oldfiles)
print,n_elements(g)
print,strlen(g)
if strlen(g) ne 0 then begin $
   print,"reading old background"& $
   b=read_gdf(bgfile) & $
   help,b & $
   b=reform(b[x0:x1,y0:y1], nx, ny, /overwrite) & $
   print,"reading old bg" & $
endif else begin $
   ;; print,(floor(n0/nmed) - 1)*nmed & $
   ;; print,n0-1 & $
   for n = 0, nmed-1 do begin $
      ;;(floor(n0/nmed) - 1)*nmed
      print, n & $
      framen = read_gdf(f[n+1]) & $
      tvscl,framen & $
      buf[n mod nmed, *, *] = (framen)[x0:x1, y0:y1] & $
   endfor & $    
   b = median(buf, /double, dimension = 1) >1 & $
   a = buf[0,*,*] & $
   c = a/(b>1) & $
   write_gdf, b, 'background.gdf' & $;;FIXME don't want it to save to real BG.
endelse  
buf = rebin(reform(b,1,nx,ny),nmed,nx,ny)
tstart = systime(1)
print,"got background"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

 

;;;;;;;;;;;;;;;;;;;;;;;;Main Loop, analyzing holograms;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for n=startfrm,endframe do begin 
   print, 'Frame', n 
   a = (read_gdf(f[n]))[5:644,*] ;;;NOTE!!!!!!!!!!!!!!!!!!
                                    ;I put the dimensions in mp_ripper
                                    ;to [656,480]. therefore I had to 
                                ; cut off the edges to get to the 
                                ; correct [640,480]. The background
                                ; is already corrected.
   aa = reform(a[x0:x1,y0:y1], nx, ny, /overwrite) 
   buf[n mod nmed, *, *] = aa 
   
   
    
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; update background ;;DON'T UPDATE BG if trapping
   ;; if n mod nmed eq 0 and n ne 0 then begin & $
   ;;    b = median(buf, /double, dimension = 1) > 1 & $
   ;;    write_gdf,b,"background.gdf"  & $
   ;; endif & $
       
   ;; if n lt startfrm then begin $;;;;;;;Pick up from where it left off
   ;;     n+=1 & $
   ;;     continue & $
   ;; endif  & $

   aa = double(aa-dc)/(b-dc>1)                              ;Normalize Image

   
   
   ;tvscl,aa                                          ;Display Image
   
   ;; wait,.5 & $
   ;; circimage = circletransform(aa,noise=.02,smooth=10,deinterlace=1) & $
   ;; ff = fastfeature(circimage,50)
   ;; print,"features:",ff
   ;; plotimage,bytscl(circimage),/iso & $
   ;; wait,1  
   ;;;;;;;;;;;;;;;;;; analyze even field;;;;;;;;;;;;;;;;;;;;;;;;
   print,"Analyzing even field 1" 
   ;help,ap,np,nm,lambda,zp,thisp,rad,mpp
   thisp = dhmfeature5(aa, ap = ap, np = np, nm = nm, lambda = lambda, $
                      threshold = featurethreshold,zguess=zp, $
                      noise=.02,smoothfactor=10,prep=thisp,$
                      rad = rad, gpu = gpu, /fit, fixalpha=fixalpha,$
                      fixzo=1,zo=0,$
                      deinterlace = 0, mpp=mpp,fixdelta=fixdelta,delta=0,$
                      fixap=fixap,fixnp=fixnp,$
                      pickn=pickn,randompos=randompos,$
                      calcfraction=calcfraction) 
   ;if doshow then plotimage, bytscl(deinterlace(aa, /odd)), /iso 
   ;;Format output
   if n_elements(thisp) gt 1 then begin 
     npts = n_elements(thisp)/18 
     thisp = reform(thisp,2,9,npts) 
     zp = thisp[0,2,0]
     ;adjust radius baised on z height
     rad = (thisp[0,2,0]*tan(35*!pi/180.)>30)<180
     ;print, thisp & $
     ;if doshow then plots, thisp[0, 0, *], thisp[0, 1, *], psym = circ() 
     p = [[[p]], [[thisp]]] 
     t = [[t], [[n, 0] # replicate(1., npts)]]
    
   endif 
   ;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;;;;;;;;;;;;;;;; analyze odd field;;;;;;;;;;;;;;;;;;;;;;;;
   print,"Analyzing odd field" 
   thisp = dhmfeature5(aa, ap = ap, np = np, nm = nm, lambda = lambda, $
                      threshold = featurethreshold,zguess=zp, $
                      noise=.02,smoothfactor=10,/quiet,prep=thisp,$
                      rad = rad, gpu = gpu, /fit, fixalpha=fixalpha,$
                      fixzo=1,zo=0,$
                      fixap=fixap,fixnp=fixnp,$
                      deinterlace = 1, mpp=mpp,fixdelta=fixdelta,delta=0,$
                      pickn=pickn,randompos=randompos,$
                      calcfraction=calcfraction) 
   if doshow then plotimage, bytscl(deinterlace(aa, /odd)), /iso 
   ;;Format output
   if n_elements(thisp) gt 1 then begin 
     npts = n_elements(thisp)/18 
     thisp = reform(thisp,2,9,npts) 
     zp = thisp[0,2,0]
     ;adjust radius baised on z height
     rad = (thisp[0,2,0]*tan(35*!pi/180.)>30)<180
     print,rad 
     ;print, thisp & $
     if doshow then plots, thisp[0, 0, *], thisp[0, 1, *], psym = circ() 
     p = [[[p]], [[thisp]]] 
     t = [[t], [[0, n] # replicate(1., npts)]]
    
   endif 
   ;;;;;;;;;;;;;;;;;;;;;;;;

   ;;;;;;;;;;;;;;;;;;;Write Results up to this frame;;;;;;;;;;;;;;;;;;;;;;
   w = where(finite(p[1, 0, *])) ; eliminate failed fits
   p = p[*, *, w] 
   t = t[*, w] 
  
  write_gdf, p, fout+"p.gdf" 
  write_gdf, t, fout+"t.gdf" 
  
  
  
  print, 'Total time:', systime(1) - tstart 
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




;;;;;;;;;;;;;;;;;;;;;Output Final Results;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
p = [[[transpose(p, [0, 2, 1])]], [[t]]] ; append time info
p = transpose(p, [0, 2, 1])
p = p[*, *, 1:*] ; clean up
write_gdf, p, fout+"full.gdf"
;end
;;;;;;;;;;;;;;;;;;;;****************end************;;;;;;;;;;;;;;;;;;;
return,p
end
