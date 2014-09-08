;+
; NAME:
;      fitspheredhm
;
; PURPOSE:
;      Measure the radius, refractive index, and three-dimensional
;      position of a colloidal sphere immersed in a dielectric 
;      medium by fitting its digital holographic microscopy (DHM)
;      image to Mie scattering theory.
;
; CATEGORY:
;      Holographic microscopy
;
; CALLING SEQUENCE:
;      params = fitspheredhm(a,p)
;
; INPUTS:
;      a : two-dimensional real-valued DHM image of sphere.
;
;      p : initial guess for fitting parameters.
;      p[0] : xp : x-coordinate of sphere's center [pixels]
;      p[1] : yp : y-coordinate [pixels]
;      p[2] : zp : z-coordinate [pixels]
;      p[3] : ap : sphere's radius [micrometers]
;      p[4] : np : sphere's refractive index
;      p[5] : amplitude : arbitrary.  1 is a reasonable value
;      p[6] : nm : refractive index of medium
;
;      Optional:
;      p[7] : kp : sphere's extinction coefficient
;      p[8] : km : medium's extinction coefficient
;
; KEYWORD PARAMETERS:
;      lambda: vacuum wavelength of illumination [micrometers].
;           Default: 0.632816
;      mpp: Length-scale calibration factor [micrometers/pixel].
;           Default: 0.135
;      precision: Accuracy with which scattering coefficients are
;           summed.  PRECISION=0.001 noticably speeds up 
;           CPU-based fits at the expense of some precision
;           in the fitting parameters.  
;           Default: PRECISION=0, retain all terms.
;
;      aplimits: [minap, maxap] limits on ap [micrometers]
;           Default: [0.05, 10.]
;
;      nplimits: [minnp, maxnp] limits on np
;           Default: [1.01 nm, 3.]
;
; KEYWORD FLAGS:
;      deinterlace: Only fit to odd (DEINTERLACE = 1)
;                   or even (DEINTERLACE = 2) scan lines.  This is
;                   useful for fitting to images that were acquired
;                   with an interlaced camera.
;
;      fixap: If set, do not allow ap to vary.
;      fixnp: If set, do not allow np or kp to vary.
;      fixnm: If set, do not allow nm or km to vary.
;      fixzp: If set, do not allow zp to vary.
;      fixalpha: If set, do not allow alpha to vary.
;
;      gpu: If set, use GPU acceleration to calculate fields on
;           systems with GPUlib installed.  Requires NVIDIA graphics
;           accelerator with CUDA support.
;
;      quiet: If set, do not show results of intermediate calculations.
;
; OUTPUTS:
;      params: Least-squares fits for the values estimated in P.
;              params[0,*]: Fit values.
;              params[1,*]: Error estimates.
;              NOTE: errors are set to 0 for parameters held constant
;              with the FIX keywords.
;
; RESTRICTIONS:
;      Becomes slower and more sensitive to accuracy of initial
;      guesses as spheres become larger.
;
; PROCEDURE:
;      Uses MPFIT by Craig Marquardt to minimize the difference
;      between the measured DHM image and the image computed by
;      SPHEREDHM.
;
; REFERENCE:
;      S. Lee, Y. Roichman, G. Yi, S. Kim, S. Yang, A. van Blaaderen,
;      P. van Oostrum and D. G. Grier,
;      Chararacterizing and tracking
;      single colloidal particles with video holographic microscopy,
;      Optics Express 15, 18275-18282 (2007)
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 4/2007.
; 05/22/2007: DGG. Added LAMBDA keyword.
; 05/26/2007: DGG. Revised to use Bohren and Huffman version of
;   SPHEREFIELD.
; 06/10/2007: DGG. Updated for more accurate BH code.
; 09/11/2007: DGG. Made nm a fitting parameter and removed NM keyword.  
;   Replaced FIXINDEX keword with FIXNM.
;   Added FIXNP keyword.  
; 11/03/2007: DGG. Changed FIXRADIUS to FIXAP.  Added FIXZP and FIXALPHA.
; 02/08/2008:  DGG. Treat coordinates as one-dimensional arrays internally
;   to eliminate repeated calls to REFORM.
;   Adopt updated syntax for SPHEREFIELD: separate x, y and z coordinates.
;   Y coordinates were incorrectly cast to float rather than double.
; 02/10/2008: DGG. Added DEINTERLACE. Small documentation fixes.
; 04/16/2008: DGG. Added MPP keyword.  Small documentation fixes.
; 10/13/2008: DGG. Added PRECISION and GPU keywords to make use of new
;   capabilities in SPHEREFIELD.
; 10/17/2008: DGG. Added LUT keyword to accelerate CPU-based fits.  
;   This required setting .STEP = 0.0001 pixel restrictions on
;   the x and y centroids in PARINFO.
; 01/15/2009: DGG. Documentation clean-ups.
; 02/14/2009: DGG. Added APLIMITS and NPLIMITS keywords.
; 03/17/2009: DGG. Added support for complex refractive indexes by
;    accounting for the particle and medium extinction coefficients,
;    kp and km.
; 03/26/2009: Fook Chiong Cheong, NYU: np and nm should be cast as
;    dcomplex rather than complex when kp or km are non-zero.
; 06/18/2010: DGG. Added COMPILE_OPT.
; 10/22/2013: DBR Added capability to scale image to account for
; interface
; 10/31/2013: DBR instead of scaling image now it adds a phase
;
; Copyright (c) 2007-2010 David G. Grier and Fook Chiong Cheong
;
; UPDATES:
;    The most recent version of this program may be obtained from
;    http://physics.nyu.edu/grierlab/software.html
; 
; LICENSE:
;    This program is free software; you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation; either version 2 of the
;    License, or (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program; if not, write to the Free Software
;    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
;    02111-1307 USA
;
;    If the Internet and WWW are still functional when you are using
;    this, you shold be able to access the GPL here: 
;    http://www.gnu.org/copyleft/gpl.html
;-

function spheredhm_f, x, y, p, $
                      lambda=lambda, mpp=mpp, $
                      precision=precision, gpu=gpu, lut=lut

COMPILE_OPT IDL2, HIDDEN

; p[0] : xp
; p[1] : yp
; p[2] : zp
; p[3] : ap
; p[4] : np
; p[5] : amplitude
; p[6] : nm
; p[7] : zo
;
; Optional:
; p[8] : kp
; p[9] : km

xx = x - p[0]
yy = y - p[1]
zp = p[2]
ap = p[3]
np = p[4]
alpha = p[5]
nm = p[6]
zo = p[7]
delta = p[8]

if n_elements(p) eq 11 then begin
   np = dcomplex(np, p[9])
   nm = dcomplex(nm, p[10])
endif

if keyword_set(lut) then begin
    rho = sqrt(xx^2 + yy^2)
    xx = findgen(1,fix(max(rho))+1)
    yy = 0. * xx
endif

;This finds the phase correction to account for the interface
rho = sqrt(xx^2 + yy^2)
ng = 1.515
phicorrect = calcphasecorrection(rho,zp,zo,nm,ng,lambda=lambda,mpp=mpp)


field = spherefield(xx, yy, zp, ap, np, nm=nm, $
                    lambda=lambda, mpp=mpp, k=k, /cartesian, $
                    precision=precision, gpu=gpu)

; interference between light scattered by the particle
; and a plane wave polarized along x and propagating along z
dhm = 1.d + 2.d * alpha * field[0,*]*exp(dcomplex(0,phicorrect+delta)) * $
                                               exp(dcomplex(0,-k*zp))

;help,phicorrect
;help,field[0,*]
;help,dhm

dhm += alpha^2 * total(field * conj(field), 1)
dhm = real_part(dhm)

if keyword_set(lut) then begin
  dhm = interpolate(dhm, rho, cubic=-0.5)
endif

return, dhm
end
 
function fitspheredhm5, a, $               ; image
                       p0, $              ; starting estimates for parameters
                       aplimits = aplimits, $ ; limits on ap [micrometers]
                       nplimits = nplimits, $ ; limits on np
                       lambda = lambda, $ ; wavelength of light [micrometers]
                       mpp = mpp, $       ; micrometers per pixel
                       fixnp = fixnp, $   ; fix particle refractive index
                       fixnm = fixnm, $   ; fix medium refractive index
                       fixap = fixap, $   ; fix particle radius
                       fixzp = fixzp, $   ; fix particle axial position
                       fixzo = fixzo, $   ; fix the height of interface
                       fixdelta = fixdelta,$ ; fix the phase offset
                       fixalpha = fixalpha, $ ; fix illumination
                       deinterlace = deinterlace, $
                       precision = precision, $ ; precision of coefficients
                       gpu = gpu, $             ; use gpu acceleration
                       lut=lut, $               ; use look-up table acceleration
                       quiet = quiet,$            ; don't print diagnostics
                       randompos = randompos,$    ;evaluate are random positions
                       calcfraction = calcfraction,$ ;fraction of points to use 
                       limitradius=limitradius               ;randompos

COMPILE_OPT IDL2

sz = size(a,/dimensions)
nx = sz[0]
ny = sz[1]
npts = nx*ny

                                ; Fixed 10_9_13 by David Ruffner using
                                ; error estimate suggested by David Grier
noise = median(abs(a-median(a,3)))
xp = p0[0]
yp = p0[1]
azia = aziavg(a,center=[xp+nx/2.,yp+ny/2.],deviates=dev)
err = abs(dev) > noise
azierr = aziavg(err,center=[xp+nx/2.,yp+ny/2.])
;plot,azia,psym=1
;oplot,azia+azierr,linestyle=1
;oplot,azia-azierr,linestyle=1
;wait,3
;err = abs(dev) > noise
;err = replicate(5.,npts) 

x = dindgen(npts) mod nx        ; coordinates of pixels
y = double(floor(dindgen(npts) / nx))
aa   = double(reform(a, npts))

if   n_elements(deinterlace) gt 0 then begin
   w   = where((y mod 2) eq (deinterlace mod 2), npts)
   x = x[w]
   y    = y[w]
   aa  = aa[w]
endif
x -= double(nx) / 2.d
y   -= double(ny) / 2.d



if n_elements(randompos) ne 0 then begin
   if n_elements(limitradius) eq 2 then begin
      print,"limiting radius"
      rho = sqrt(x^2+y^2)
      b = where(rho lt limitradius[1] and rho gt limitradius[0])
      x = x[b]
      y = y[b]
      aa=aa[b]
      npts = n_elements(x)
   endif
   if n_elements(calcfraction) ne 0 then begin
      numpts = round(npts*calcfraction)
   endif
   indices = Round( RandomU(seed, numpts) * npts )
   x = x[indices]
   y = y[indices]
   aa = aa[indices]
   if not keyword_set(quiet) then $
       plotimage, bytscl(a), /iso & $
       oplot, x+double(nx) / 2.d, y+double(ny) / 2.d, psym = circ(/fill) & $
endif

if n_elements(lambda) ne 1 then $
   lambda = 0.632816d           ; HeNe wavelength (micrometers)

if n_elements(mpp) ne 1 then $
   mpp = 0.101                  ; Zeiss rig

if n_elements(precision) ne 1 then $
   precision = 0.               ; Keep all scattering coefficients

if n_elements(gpu) ne 1 then $
   gpu = 0.

if n_elements(lut) ne 1 then $
   lut = 0

nparams = n_elements(p0)
parinfo = replicate({limited:[0,0], limits:[0.,0.], fixed:0, step:0.}, nparams)
;; Restrictions on fitting parameters
; xp and yp: overly small steps prevent convergence of LUT-based fits.
parinfo[0:2].step = 0.0001      ; No apparent harm for regular or GPU fits.
; zp: No restrictions
; ap: Radius must be positive
parinfo[3].limited[0] = 1
parinfo[3].limits[0] = 0.05
parinfo[3].limited[1] = 1
parinfo[3].limits[1] = 10.
if n_elements(aplimits) eq 2 then $
   parinfo[3].limits = aplimits
; np: Refractive index of particle
parinfo[4].limited[0] = 1
parinfo[4].limits[0] = 0.01*p0[6] ; FIXME what about low-index particles?
parinfo[4].limited[1] = 1
parinfo[4].limits[1] = 3.0      ; a bit more than titania
if n_elements(nplimits) eq 2 then $
   parinfo[4].limits = nplimits
; alpha: Illumination at particle
parinfo[5].limited[0] = 1
parinfo[5].limits[0] = 0.       ; cannot be negative
; nm: Refractive index of medium: No restrictions
; Flags to prevent adjusting parameter values
parinfo[2].fixed = keyword_set(fixzp)
parinfo[3].fixed = keyword_set(fixap)
parinfo[4].fixed = keyword_set(fixnp)
parinfo[5].fixed = keyword_set(fixalpha)
parinfo[6].fixed = keyword_set(fixnm)

parinfo[7].fixed = keyword_set(fixzo)
parinfo[7].limited[*] = 1
parinfo[7].limits[0] = 0
parinfo[7].limits[1] = p0[2]*1.515/p0[6]
parinfo[8].fixed = keyword_set(fixdelta)
parinfo[8].limited[*] = 1
parinfo[8].limits[0] = -!pi
parinfo[8].limits[1] = !pi

if nparams eq 11 then begin
   parinfo[9].fixed = keyword_set(fixnp)
   parinfo[10].fixed = keyword_set(fixnm)
   parinfo[9:10].limited[*] = 1
   parinfo[9:10].limits[0] = 0.
   parinfo[9:10].limits[1] = 5.
endif

; parameters passed to the fitting function
argv = {lambda:lambda, mpp:mpp, precision:precision, gpu:gpu, lut:lut}

; errors from fit
perror = fltarr(nparams)

; perform fit
p = mpfit2dfun("spheredhm_f", x, y, aa, err, p0, functargs=argv, $
               parinfo = parinfo, /fastnorm, $
               perror=perror, bestnorm=bestnorm, dof=dof, $
               status = status, errmsg = errmsg, quiet=quiet)




if status le 0 then begin 
   message, errmsg, /inf
   return, -1
endif

; failure?
if n_elements(p) eq 1 then begin
   message, "MPFIT2DFUN did not return a result",  /inf
   return, -1
endif

; success
; rescale fit uncertainties into error estimates
dp = perror* sqrt(bestnorm/dof)

return, [transpose(p),transpose(dp)]
end
