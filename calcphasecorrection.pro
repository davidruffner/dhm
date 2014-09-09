;+
; NAME:
; calcphasecorrection
;
; PURPOSE:
; calculate the correction to the phase of the scattered field of a
; particle in digital holographic microscopy(DHM) due to the interface
; between the water and the glass.
;
; CATEGORY:
; holographic video microscopy
;
; CALLING SEQUENCE:
; phase =
; calcphasecorrection(rgs,zdhm,zo,nw,ng,lambda=lambda,mpp=mpp)
;
; INPUTS:
; rgs: radius of point away from the particle in the focal plane,
; in pix
; zdhm: height of the particle estimated by DHM, in pix
; zo: distance from the focal plane to the oil/glass-water, in pix
; interface
; nw: refractive index of the water
; ng: refractive index of the oil/glass
;
; KEYWORD PARAMETERS:
; lambda: vacuum wavelength of the imaging light
; mpp: microns per pixel
;
; OUTPUT:
; phase: array containing correcting phase for each pixel
;
; PROCEDURE:
; The small angle approximation is used to find where a ray ending
; at a point in the focal plane would cross the interface. Then
; using this point the phase is calculated.
;
; MODIFICATION HISTORY:
; 11/08/2013: Written by David Ruffner, New York University
;-
function calcphasecorrection,rgs,zdhm,zo,nw,ng,lambda=lambda,mpp=mpp
if n_elements(lambda) eq 0 then lambda=.447
if n_elements(mpp) eq 0 then mpp = .135
nw = real_part(nw)
ng = real_part(ng)
zg = zo ;the oil and the glass have the same refractive index
zactual = zdhm+zg*(1-nw/ng) ; The height we estimate the particle to actually be
zw = zactual-zg
npts = n_elements(rgs)
if zw lt 0 then begin
print,"zo is too large!"
return,fltarr(npts)
endif

;Switch to micrometers
zactualum = zactual*mpp
zwum = zw*mpp
zgum = zg*mpp
rgsum = rgs*mpp
zdhmum= zdhm*mpp

;Calculate where the ray would intersect the interface
alpha = 1+nw*zg/(ng*zw)
rwsum = rgsum/alpha

;Calculate the phase of ray at focal plane
lwsum = sqrt(zwum^2+rwsum^2) ;We use square roots here which makes it
; more accurate
lgsum = Sqrt(zgum^2 + (rwsum - rgsum)^2)

phiscat = (2*!pi/lambda)*(nw*lwsum+ng*lgsum)
phisph = (2*!pi/lambda)*((nw*sqrt(rgsum^2 + zdhmum^2))$
+(zgum*ng - zgum*nw^2/ng))

phicorrect = phiscat-phisph

return,reform(phicorrect)
end
