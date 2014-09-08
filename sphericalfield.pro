;+
; NAME:
;       sphericalfield
;
; PURPOSE:
;       Calculates the complex electric field defined by an array of scattering
;       coefficients.
;
; CATEGORY:
;       Holography, light scattering, microscopy
;
; CALLING SEQUENCE:
;       field = sphericalfield(x, y, z, a, lambda, $
;                           mpp = mpp)
;
; INPUTS:
;       x: [npts] array of pixel coordinates [pixels]
;       y: [npts] array of pixel coordinates [pixels]
;       z: If field is required in a single plane, then
;          z is the plane's distance from the sphere's center
;          [pixels].
;          Otherwise, z is an [npts] array of coordinates.
;
;       NOTE: Ideally, x, y and z should be double precision.
;             This is left to the calling program for efficiency.
;
;       a: [2,nc] array of a and b scattering coefficients, where
;          nc is the number of terms required for convergence.
;
;       lambda: wavelength of light in medium [pixels]
;
; KEYWORD FLAGS:
;       cartesian: If set, return field components in Cartesian
;           coordinates.  Default: Spherical polar coordinates
;
; OUTPUTS:
;       field: [3,npts] complex values of field at the positions r.
;              [0,*]: r component
;              [1,*]: theta component
;              [2,*]: phi component
;
;              If CARTESIAN is set:
;              [0,*]: x component (incident polarization)
;              [1,*]: y component (transverse component)
;              [2,*]: z component (axial component, relative to
;              incident beam).
;
; REFERENCE:
;   1. Adapted from Chapter 4 in
;      C. F. Bohren and D. R. Huffman, 
;      Absorption and Scattering of Light by Small Particles,
;      (New York, Wiley, 1983).
;   2. W. J. Wiscombe,
;      Improved Mie scattering algorithms,
;      Applied Optics 19, 1505-1509 (1980).
;
; MODIFICATION HISTORY:
; Written by David G. Grier, New York University, 5/2007
; 6/9/2007: DGG finally read Section 4.8 in Bohren and Huffman about
;    numerical stability of the recursions used to compute the scattering
;    coefficients.  Feh.  Result is a total rewrite.
; 6/20/2007: DGG Calculate \tau_n(\cos\theta) and \pi_n(\cos\theta)
;    according to recurrence relations in 
;    W. J. Wiscombe, Appl. Opt. 19, 1505-1509 (1980).
;    This is supposed to improve numerical accuracy.
; 2/8/2008: DGG. Replaced single [3,npts] array of input coordinates
;    with two [npts] arrays for x and y, and a separate input for z.
;    Eliminated double() call for coordinates.  Z may have 1 element or
;    npts elements. Small documentation fixes.
; 4/3/2008: Bo Sun (Sephiroth), NYU: Calculate Lorenz-Mie a and b
;    coefficients using continued fractions rather than recursion.
;    Osman Akcakir from Arryx pointed out that the results are
;    more accurate in extreme cases.  Method described in
;    William J. Lentz, "Generating Bessel functions in Mie scattering
;    calculations using continued fractions," Appl. Opt. 15, 668-671
;    (1976).
; 4/4/2008: DGG small code clean-ups and documentation.  Added
;    RECURSIVE keyword for backward compatibility in computing a and b
;    coefficients.
; 4/11/2008: Sephiroth: Corrected small error in jump code for
;    repeated fractions in Mie coefficients.
; 6/25/2008: DGG Don't clobber x coordinate input values.
; 10/9/2008: DGG adapted from SPHEREFIELD by separating out
;    calculation of scattering coefficients, a_n and b_n.  This
;    is therefore more general, and can be replaced more
;    readily with a GPU-accelerated version.
; 10/13/2008: DGG eliminated RECURSIVE keyword.
;
; Copyright (c) 2007-2010 Bo Sun and David G. Grier
;-

function sphericalfield, x_, y_, z_, ab, lambda, $
                         cartesian = cartesian ; project to cartesian coordinates

npts = n_elements(x_)
nc = n_elements(ab[0,*])-1      ; number of terms required for convergence

k = 2.d * !dpi / lambda         ; wavenumber in medium [pixel^-1]

ci = dcomplex(0,1)

; convert to spherical coordinates centered on the sphere.
; (r, theta, phi) is the spherical coordinate of the pixel
; at (x,y) in the imaging plane at distance z from the
; center of the sphere.
rho   = sqrt(x_^2 + y_^2)
r     = sqrt(rho^2 + z_^2)
theta = atan(rho, z_)
phi   = atan(y_, x_)
costheta = cos(theta)
sintheta = sin(theta)
cosphi = cos(phi)
sinphi = sin(phi)

kr = k*r                        ; reduced radial coordinate

; starting points for recursive function evaluation ...
; ... Riccati-Bessel radial functions, page 478
sinkr = sin(kr)
coskr = cos(kr)
xi_nm2 = dcomplex(coskr, sinkr) ; \xi_{-1}(kr)
xi_nm1 = dcomplex(sinkr,-coskr) ; \xi_0(kr)

; ... angular functions (4.47), page 95
pi_nm1 = 0.d                    ; \pi_0(\cos\theta)
pi_n   = 1.d                    ; \pi_1(\cos\theta)

; storage for vector spherical harmonics: [r,theta,phi]
Mo1n = dcomplexarr(3,npts)
Ne1n = dcomplexarr(3,npts)

; storage for scattered field
Es = dcomplexarr(3,npts)

; Compute field by summing multipole contributions
for n = 1.d, nc do begin

; upward recurrences ...
; ... Legendre factor (4.47)
; Method described by Wiscombe (1980)
    swisc = pi_n * costheta 
    twisc = swisc - pi_nm1
    tau_n = n * twisc - pi_nm1  ; \tau_n(\cos\theta)

; ... Riccati-Bessel function, page 478
    xi_n   = (2.d*n - 1.d) * xi_nm1 / kr - xi_nm2    ; \xi_n(kr)

; vector spherical harmonics (4.50)
;   Mo1n[0,*] = 0.d             ; no radial component
    Mo1n[1,*] = pi_n * xi_n     ; ... divided by cosphi/kr
    Mo1n[2,*] = -tau_n * xi_n   ; ... divided by sinphi/kr

    dn = (n * xi_n)/kr - xi_nm1
    Ne1n[0,*] = n*(n + 1.d) * pi_n * xi_n ; ... divided by cosphi sintheta/kr^2
    Ne1n[1,*] = -tau_n * dn     ; ... divided by cosphi/kr
    Ne1n[2,*] = pi_n  * dn      ; ... divided by sinphi/kr

; prefactor, page 93
    En = ci^n * (2.d*n + 1.d)/ n / (n + 1.d)

; the scattered field in spherical coordinates (4.45)
    Es += En * (ci * ab[0,n] * Ne1n - ab[1,n] * Mo1n)

; upward recurrences ...
; ... angular functions (4.47)
; Method described by Wiscombe (1980)
    pi_nm1 = pi_n
    pi_n = swisc + (n + 1.d) * twisc / n

; ... Riccati-Bessel function
    xi_nm2 = xi_nm1
    xi_nm1 = xi_n
endfor

; geometric factors were divided out of the vector
; spherical harmonics for accuracy and efficiency ...
; ... put them back at the end.
Es[0,*] *= cosphi * sintheta / kr^2
Es[1,*] *= cosphi / kr
Es[2,*] *= sinphi / kr

; By default, the scattered wave is returned in spherical
; coordinates.  Project components onto Cartesian coordinates.
; Assumes that the incident wave propagates along z and 
; is linearly polarized along x
if keyword_set(cartesian) then begin
    Ec = Es
    Ec[0,*] =  Es[0,*] * sintheta * cosphi
    Ec[0,*] += Es[1,*] * costheta * cosphi
    Ec[0,*] -= Es[2,*] * sinphi

    Ec[1,*] =  Es[0,*] * sintheta * sinphi
    Ec[1,*] += Es[1,*] * costheta * sinphi
    Ec[1,*] += Es[2,*] * cosphi

    Ec[2,*] =  Es[0,*] * costheta - Es[1,*] * sintheta

    return, Ec
endif

return, Es
end
