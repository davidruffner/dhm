;+
; NAME:
;    dhmtool
;
; PURPOSE:
;    Given a holographic video microscope snapshot containing
;    the image of at least one sphere, perform a Lorenz-Mie fit 
;    the region around the sphere.  Return the refined fits for 
;    the sphere's radius, refractive index and position relative 
;    to the image's origin.
;
; CATEGORY:
;    Holographic video microscopy
;
; CALLING SEQUENCE:
;    p = dhmtool(a, zp, ap, np, alpha, nm)
;
; INPUTS:
;    a: holographic snapshot containing images of spheres
;    zp: sphere's axial displacement from the focal plane [pixels]
;    ap: sphere's radius [micrometers]
;    np: sphere's refractive index.
;    alpha: estimate for normalized illumination at sphere.
;        1 is a good starting point.
;    nm: refractive index of medium.
;
; KEYWORD PARAMETERS:
;    rp: [xp, yp] sphere's centroid [pixels].  
;        Default: Interactive, if not set.
;
;    rad: range of interest [pixels].
;        Default: 50.
;
;    lambda: vacuum wavelength of light [micrometers]
;        Default: 0.6328
;
;    mpp: micrometers per pixel.
;        Default: 0.135
;
;    deinterlace: Fit only even field (DEINTERLACE = 0) or 
;        odd field (DEINTERLACE = 1) of (interlaced) hologram.
;        Default: No deinterlacing
;
;    aplimits: [minap, maxap] limits on ap [micrometers]
;           Default: [0.05, 10.]
;
;    nplimits: [minnp, maxnp] limits on np
;           Default: [1.01 nm, 3.]
;
; KEYWORD FLAGS:
;    fixnp: do not adjust np; p[1,4] = 0
;    fixalpha: do not adjust alpha; p[1,5] = 0
;    fixnm: do not adjust nm; p[1,6] = 0
;
;    lut: use look-up table acceleration on CPU-based calculations.
;    gpu: use GPU acceleration, if available.
;
;    quiet: do not provide diagnostic output.
;
; OUTPUTS:
;    p: [2,7] array of fitting parameters and error estimates
;    p[0,0]: xp [pixels]
;    p[0,1]: yp [pixels]
;    p[0,2]: zp [pixels]
;    p[0,3]: ap [micrometers]
;    p[0,4]: np
;    p[0,5]: alpha
;    p[0,6]: nm
;
;    p[1,*] contains associated error estimates.
;
; SIDE EFFECTS:
;    If RC is not provided, then DHMTOOL presents an image
;    and requires the user to click on the sphere's center.
;    Unless QUIET is set, DHMTOOL provides diagnostic output.
;
; RESTRICTIONS:
;    Requires a normalized hologram.
;
; PROCEDURE:
;    Crop region around the estimated centroid.  Perform fit
;    on cropped image using FITSPHEREDHM.  Adjust fit parameters 
;    to account for cropping.
;
; MODIFICATION HISTORY:
; 01/15/2009 Written by David G. Grier, New York University.
; This formalizes a routine that has been in use since 2007.
; 02/14/2009 DGG. Added APLIMITS and NPLIMITS keywords.
;    Documentation fixes.
;
; Copyright (c) 2007-2010 David G. Grier.
;- 

function dhmtool5, a, $               ; hologram
                  zp, $              ; z offset of particle [pixel]
                  ap, $              ; radius of particle [micrometers]
                  np, $              ; refractive index of particle
                  alpha, $           ; fraction of scattered light
                  nm, $              ; refractive index of medium
                  rc = rc, $         ; centroid [pixels]
                  rad = rad, $       ; range of interest [pixels]
                  lambda = lambda, $ ; wavelength of light [micrometers]
                  mpp = mpp, $       ; micrometers per pixel
                  aplimits = aplimits, $
                  nplimits = nplimits, $
                  fixnp = fixnp, $
                  fixnm = fixnm, $
                  fixalpha = fixalpha, $
                  fixdelta = fixdelta, $
                  fixap = fixap, $
                  fixzo = fixzo, $
                  deinterlace = deinter, $
                  lut = lut, gpu = gpu, $ ; acceleration
                  quiet = quiet,$           ; minimize output
                  randompos = randompos,$    ;evaluate are random positions
                  calcfraction = calcfraction,$ ;fraction of points to use 
                                                   ;randompos 
                  limitradius=limitradius, $
                  zo = zo,$                  ;distance from the focal plane to the                                          ;glass water interface
                  delta=delta               ;phase offset

if n_elements(rad) lt 1 then rad = 50
if n_elements(zo) eq 0 then zo=0
if n_elements(delta) eq 0 then delta=0

; estimate center of scattering pattern, if not provided
if n_elements(rc) ne 2 then begin
    plotimage, bytscl(a), /iso
    print, "Click on center of sphere"
    cursor, x, y, /data
    xc = x
    yc = y
    rc = [xc, yc]
endif else begin
    xc = rc[0]
    yc = rc[1]
endelse

; clip hologram to region around center
sz = size(a, /dimensions)
x0 = round(xc - rad) > 0
x1 = round(xc + rad) < sz[0]-1
y0 = round(yc - rad) > 0
y1 = round(yc + rad) < sz[1]-1

ac = a[x0:x1, y0:y1]            ; cropped image

chatty = not keyword_set(quiet)

;; if chatty then $
;;    if n_elements(deinter) ge 1 then $
;;       plotimage, bytscl(deinterlace(ac, odd=deinter)), /iso $
;;    else $
;;       plotimage, bytscl(ac), /iso

; dimensions of clipped image
nx = x1 - x0 + 1
ny = y1 - y0 + 1

; particle position relative to center of clipped image
xp = xc - x0 - nx/2.
yp = yc - y0 - ny/2.
;; print,"rc",rc
;; print,"x0",x0
;; print,"nx/2",nx/2
; show initial estimates
if chatty then begin
   if n_elements(deinter) ge 1 then $
      tvscl, deinterlace(ac, odd=deinter) $ 
   else $
      tvscl, ac
   b = spheredhm5([xp,yp,zp], ap, np, nm, [nx, ny], $
                 alpha = alpha, lambda = lambda, mpp = mpp,zo=zo)
   tvscl, b, nx+1, 0
endif 

;initialize guesses
pin = [xp, yp, zp, ap, np, alpha, nm,zo,delta]
print,"pin",pin

; perform fit

if n_elements(deinter) ge 1 then $
   p = fitspheredhm5(ac, pin, $
                    lambda = lambda, mpp = mpp, $
                    aplimits = aplimits, nplimits = nplimits, $
                    fixnm = fixnm, fixnp = fixnp, $
                    fixap = fixap, fixalpha = fixalpha,fixzo=fixzo, $
                    fixdelta=fixdelta,$
                    lut = lut, gpu = gpu, $
                    deinterlace = y0 + deinter, quiet=quiet,$
                    randompos = randompos,calcfraction = calcfraction,$
                    limitradius=limitradius) $
else $
   p = fitspheredhm5(ac, pin, $
                    lambda = lambda, mpp = mpp, $
                    aplimits = aplimits, nplimits = nplimits, $
                    fixnm = fixnm, fixnp = fixnp, $
                    fixap = fixap, fixalpha = fixalpha,fixzo=fixzo, $
                    fixdelta=fixdelta,$
                    lut = lut, gpu = gpu, quiet=quiet,randompos = randompos,$
                    calcfraction = calcfraction,limitradius=limitradius)

if n_elements(p) eq 1 then return, -1

if chatty then begin
   noise = median(abs(ac-median(a,3)))
   xp = p[0,0]
   yp = p[0,1]
   azia = aziavg(ac,center=[xp+nx/2.,yp+ny/2.],deviates=dev)
   err = abs(dev) > noise
   azierr = aziavg(err,center=[xp+nx/2.,yp+ny/2.])
   npts = n_elements(azia)
   rho = findgen(npts)
   
   plot,rho,azia,psym=1
   oplot,rho,azia+azierr,linestyle=1
   oplot,rho,azia-azierr,linestyle=1
   b = spheredhm5(p[0,0:2], p[0,3], p[0,4], nm, [nx, ny], $
                   alpha = p[0,5], lambda = lambda, mpp = mpp,zo=p[0,7])
   fit = aziavg(b,center=[xp+nx/2.,yp+ny/2.]) 
   print,"plotting in dhmtool"
   ;; fit = spheredhmprofile5(rho,p[0,2],p[0,3],$
   ;;                         p[0,4],nm,p[0,5],$
   ;;                        lambda=lambda,mpp=mpp,zo=p[0,7],delta=p[0,8])
   oplot,rho,fit
   
   tvscl, [ac,b,ac-b+1.]
   
   print, "Height: "+strtrim(p[0,2],2)+" +/- "+strtrim(p[1,2],2)+" [pixels]"
   print, "Radius: "+strtrim(p[0,3],2)+" +/- "+strtrim(p[1,3],2)+" [micron]"
   print, "Index : "+strtrim(p[0,4],2)+" +/- "+strtrim(p[1,4],2)
endif
;write_gdf,p,"rawp.gdf"
p[0,0:1] += [x0+nx/2., y0+ny/2.] 

return, p
end
