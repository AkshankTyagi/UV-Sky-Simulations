; GDL Version 1.0.0-rc.3 (darwin arm64 m64)
; Journal File for jayanth@BigMac.local
; Working directory: /Users/jayanth/Dropbox/user/data/castelli
; Date: Sun Oct 29 21:11:00 2023

d=mrdfits('ckm10/ckm10_18000.fits',0,hdr)
;MRDFITS: Null image, NAXIS=0
d=mrdfits('ckm10/ckm10_18000.fits',1,hdr)
;MRDFITS: Binary table.  12 columns by  1221 rows.
help,/st,d
plot,d.wavelength,d.g25
plot,d.wavelength,d.g25,xrange=[0,1000]
plot,d.wavelength,d.g25,xrange=[0,2000]
plot,d.wavelength,d.g25,xrange=[0,3000]
d1=mrdfits('ckp00/ckp00_20000.fits',1,h1)
;MRDFITS: Binary table.  12 columns by  1221 rows.
