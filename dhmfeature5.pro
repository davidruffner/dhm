;+
; NAME:
;    dhmfeature
;
; PURPOSE:
;    Identify, locate and optionally characterize spheres in
;    normalized holographic video microscopy images.
;
; CATEGORY:
;    Image analysis, holographic video microscopy, feature detection
;
; CALLING SEQUENCE:
;    p = dhmfeature(image)
;
; INPUTS:
;    image: two-dimensional normalized holographic microscopy image.
;
; KEYWORD PARAMETERS:
;    Parameters for 2D feature location
;    noise: Estimate for the RMS additive noise at each pixel.
;        Default: 0.1
;    smoothfactor: Smoothing to reduce effects of speckle.
;        Default: 5 [pixels]
;    threshold: Minimum weighting for a pixel to be considered part
;        of a feature.
;        Default: 100
;    prep : previous values of p, used to provide a starting guess,
;           it does a simple minimum distance to associate each
;           particle in the new frame with particles in the old.
;
;    Parameters for 3D feature location and characterization:
;    rad: Radius over which to fit a given feature [pixels]
;        Default: 100 
;    ap: ballpark radius of sphere [micrometers]
;        Default: 0.75
;    np: ballpark refractive index of particle.
;        Default: 1.59 (Polystyrene)
;    nm: refractive index of medium.
;        Default: water at room temperature for the given wavelength.
;    lambda: vacuum wavelength of light [micrometers]
;        Default: 0.532
;    mpp: length calibration [micrometers per pixel]
;        Default: 0.101
;
;    limitap: [minap, maxap] limits on ap [micrometers]
;        Default: [0.05, 10.]
;
;    limitnp: [minnp, maxnp] limits on np
;        Default: [1.01 nm, 3.]
;
; KEYWORD FLAGS:
;    fit: If set, perform fits for 3D feature location and
;        characterization.  Otherwise perform 2D feature location.
;    deinterlace: If set to an even number, consider only even
;        lines in feature extraction.  Similarly, if set to an
;        odd number.
;    gpu: If set, use hardware acceleration through GPULIB.
;        Requires appropriate hardware and installation of GPULIB.
;
;    quiet: If set, report only errors
;
; OUTPUTS:
;    p: [7,2,npts] array of data for the features found in the image.
;       npts: Number of features found
;       p[0,*,n]: [xp, dxp] : x location and error [pixels]
;       p[1,*,n]: [yp, dyp] : y location and error [pixels]
;       p[2,*,n]: [zp, dzp] : z location and error [pixels]
;       p[3,*,n]: [ap, dap] : radius and error [micrometers]
;       p[4,*,n]: [np, dnp] : refractive index and error
;       p[5,*,n]: [alpha, dalpha] : relative illumination at rp, and error
;       p[6,0,n]: nm: refractive index of medium
;      
; PROCEDURE:
;     circletransform highlights positions of sphere in the plane.
;     fastfeature identifies initial particle positions in the plane.
;     fitspheredhm1d provides rough starting estimate for zp, ap and
;     np.
;     dhmtool refines these estimates.
;
; MODIFICATION HISTORY:
; 10/07/08: Written by FC Cheong and David G. Grier, New York University.
; 11/07/08: Modified by FC Cheong to include keyword FIT
; 12/07/08: Modified by FC Cheong to include Rayleigh Sommerfeld
;    approximation
; 01/27/2009: DGG: Major overhaul.  Use fitspheredhm1d rather than
;    fresneldhm to estimate zp, ap and np.  Implemented DEINTERLACE
;    First draft of algorithm to filter against "bad" features.
; 02/08/2009: DGG: Avoid fitting to badly estiamted features.
;    This prevents lock-ups in dhmtool.  Major documentation overhaul.
; 02/14/2009: DGG: Added LIMITAP and LIMITNP keywords.
; 04/21/2011: DBR: Added PREP keyword
;
; Copyright (c) 2008-2010 Fook Chiong Cheong and David G. Grier.
;-

function dhmfeature5, a, $
                     noise = noise, $
                     smoothfactor = smoothfactor, $
                     threshold = threshold, $
                     fit = fit, $
                     rad = rad, $
                     ap = ap, np = np, nm = nm, $
                     limitap = limitap, limitnp = limitnp, $
                     lambda = lambda, mpp = mpp, $
                     gpu = gpu, lut = lut, $
                     deinterlace = deinter, $
                     quiet = quiet, debug = debug, zguess = zguess, $
                     prep = prep, $
                     brightcut = brightcut,$
                     pickn=pickn,$
                     randompos = randompos,$    ;evaluate are random positions
                     calcfraction = calcfraction,$ ;fraction of points to use 
                                                   ;randompos
                     limitradius=limitradius,$
                     fixalpha=fixalpha,$
                     fixdelta=fixdelta,$
                     fixzo = fixzo,$
                     fixap = fixap,$
                     fixnp = fixnp,$
                     zo = zo,$       ;distance from the focal plane to the                                          ;glass water interface
                     delta=delta

if n_elements(noise) ne 1 then noise = 0.1
if n_elements(threshold) ne 1 then threshold = 100
if n_elements(brightcut) ne 1 then brightcut = 0

minz = 100.
maxz = 400. 

debug = keyword_set(debug)

sz = size(a)
if sz[0] ne 2 then begin
   message, "USAGE: IDL> p = dhmfeature(a)", /inf
   message, "A must be a two-dimensional normalized DHM image"
   return, -1
endif
nx = sz[1]
ny = sz[2]

;if n_elements(deinter) eq 1 then $
;   a = deinterlace(ain, odd=deinter) $
;else $
;   a = ain

;;;;;;;;Use circle transform and then fastfeature to find particles;;;;;
b = circletransform(a, noise = noise, smoothfactor = smoothfactor, $
                    deinterlace = deinter) ; find circular features
if debug then plotimage, bytscl(b), /iso

thisp = fastfeature(b, threshold,brightcut=brightcut) ; find peaks
;print,thisp
npts = n_elements(thisp) / 3



if npts le 0 then return, -1

if n_elements(pickn) ne 0 then begin
   featindex = sort(thisp[2,*])
   thisp = thisp[*,reverse(featindex)]
   thisp = thisp[*,0:pickn-1]
   npts=pickn
endif
if not keyword_set(quiet) then print, npts, " features found"
;if n_elements(zo) eq 1 then numparams = 8 else numparams = 7
numparams = 9
p = reform(dblarr(2, numparams,npts),2,numparams,npts)
p[0,0:1,*] = thisp[0:1,*]

if not keyword_set(quiet) then begin
   plotimage, bytscl(a), /iso
   oplot, p[0,0,*], p[0,1,*], psym = circ(/fill)
endif

if keyword_set(fit) then begin
   ; use fitspheredhm to fit the sphere
   if n_elements(rad) ne 1 then $
      rad = 100                 ; range over which to fit [pixels]
   if n_elements(ap) ne 1 then $
      ap = 0.75                 ; radius of sphere [micrometers]
   if n_elements(np) ne 1 then $
      np = 1.59                 ; refractive index of sphere [polystyrene]
   if n_elements(lambda) ne 1 then $
      lambda = 0.532            ; laser wavelength [micrometers]
   if n_elements(nm) ne 1 then $
      nm = refractiveindex(lambda, 24.) ; refractive index at room temperature
   if n_elements(mpp) ne 1 then $
      mpp = 0.101               ; microns per pixel

   alpha = 1.
      
   for j = 0, npts-1 do begin
      ;don't bother to fit a noisy image
      as = azistd(a, aa, center = p[0, 0:1, j], deinterlace = deinter)
      if n_elements(as) lt rad then begin
         if not keyword_set(quiet) then print,"must be on the edge"
         if n_elements(as) lt 30 then continue else rad = n_elements(as)-1
      endif
      crud = total(as[1:rad])/total((aa[1:rad]-1.)^2) 
      if not keyword_set(quiet) then print, j,"  crud",crud
      if crud gt 45. then continue
       
 
             ;;NOTE the limiting value is usually 4

      ; estimate starting parameters for fits
      aa = aziavg(a, center = p[0, 0:1, j], rad = rad, deinterlace = deinter)
      if n_elements(aa) le 1 then begin
         if not keyword_set(quiet) then print, "po=", p[0, 0:1, j]
         continue
      endif

      ;; assuming the estimates for ap and np are close enough,

          ;;Using a guess for z instead of the dhmzscan
      ;;zp = dhmzscan(aa, ap, np, nm, lambda=lambda, mpp = mpp)
      if n_elements(zguess) eq 0 then begin
         zp =  53
      endif else zp = zguess

      ;;;Estimate for z based on previous fits;;;;;
      nprep = n_elements(prep)/(2*numparams);if prep=-1 (i.e. there was a bad fit)
                                 ;it doesn't go through this process
      lmin = 0
      dr = 10000
      if nprep ne 0 then begin
         if not keyword_set(quiet) then $
              print,"using previous results to get zguess"
         ;print,"prep"
         ;print,prep
         ;print,"these x, y values"
         ;print,reform(p[0,0:1,j])
                  
         for l=0,nprep-1 do begin
            ;print,"l",l
            drtemp = Sqrt(total((p[0,0:1,j]-prep[0,0:1,l])^2.))
            if drtemp lt dr then begin
               lmin = l
               dr = drtemp
            endif
         endfor
         zp = prep[0,2,lmin]
         if not keyword_set(quiet) then print,"zp guess",zp
      endif
           
     
            
      ; otherwise fit feature
      thisp = dhmtool5(a, zp, ap, np, alpha, nm, $
                      rad=rad, rc=p[0,0:1,j], $
                      lambda = lambda, mpp = mpp, $
                      aplimits = limitap, nplimits = limitnp, $
                      /fixnm, fixalpha=fixalpha,$
                      fixdelta=fixdelta,fixzo=fixzo,$
                      gpu = gpu, lut = lut, $
                      deinterlace = deinter, $
                      quiet = quiet,randompos = randompos,$
                      calcfraction = calcfraction,limitradius=limitradius,$
                      zo=zo,delta=delta,fixap=fixap,fixnp=fixnp)
      if n_elements(thisp) eq 2*numparams then $
         p[*,*,j] = thisp
      if not keyword_set(quiet) then begin
         print, 'particle ', j
         print, "p fit",p[*,*,j]
      endif
   endfor
   w = where(finite(p[1, 0, *]) and p[0, 2, *] ne 0, ngood)
   if ngood lt npts then begin
      message, "some fits failed.  Continuing ...", /inf
      write_gdf, a, "badhologram.gdf"
      write_gdf, p, "badparams.gdf"
   endif
   if ngood le 0 then return, -1
   p = p[*, *, w]
endif

return, p
end
