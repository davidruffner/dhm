;program to rip a movie and to calculate the
;background

;Keywords
;
;fout:      The name to be used to save the background file.
;
;path:      name and path of video file to be ripped. If specified
;           it will rip the video
;
;fname:     name of folder with the ripped frames
;
;startfrm:  first frame number to begin taking background
;
;endframe:  last frame number to end background
;
;step:      frame number sampling interval for background

pro getbackground,frames=frames,fout=fout,path=path,fname = fname,$
                 startfrm=startfrm,endframe=endframe, step=step
;;Parameters
w=656
h=480
if n_elements(path) eq 0 then path = '../VTS_01_1.VOB'
if n_elements(fname) eq 0 then fname = 'vid01frames'
if n_elements(fout) eq 0 then fout="background.gdf"
help,frames
f = file_search(fname+"/*/*/*.gdf")
if n_elements(frames) ne 0 then begin
   print,"using the inputted frames" 
   f = frames
endif
if n_elements(f) eq 1 then begin $
   if f eq '' then begin $
      print,"ripping video.." & $
      mp_ripper,path,fname,dimensions = [w,h], /greyscale & $
      f = file_search(fname+"/*/*/*.gdf") & $
   endif & $
endif
 


 
nfrms = n_elements(f)

;;First get a background
nbg = 200 < nfrms 
stack = bytarr(w,h)
if n_elements(startfrm) eq 0 then startfrm = 0
if n_elements(endframe) eq 0 then endframe = startfrm+nbg
if endframe eq -1 then endframe= nfrms-1
if n_elements(step) eq 0 then step = 5
Print,"Getting BG"
print,"endframe",endframe
for i=startfrm,endframe,step do begin $
   a = read_gdf(f[i]) & $
   ;plotimage,bytscl(a),/iso & $
   ;help,a & $
   stack = [[[stack]],[[a]]] & $
   print,i & $
   plotimage,bytscl(a),/iso & $
   endfor
help,stack
stack = stack[5:644,*,1:-1]
bg = median(stack,dimension=3)
help,bg
plotimage,bytscl(bg),/iso
write_gdf,bg,fout

a = (read_gdf(f[startfrm]))[5:644,*] ;read_frame1


;write_image,"frame"+strtrim(startfrm,2)+".png","png",bytscl(a)
;write_gdf,a,"frame"+strtrim(startfrm,2)+".gdf"

;write_image,"frame"+strtrim(startfrm,2)+"clean.png","png",bytscl(a/(bg>1))
;write_gdf,a/(bg>1),"frame"+strtrim(startfrm,2)+"clean.gdf"


end


