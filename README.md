DHM

A repository of IDL routines for analyzing videos using digital holographic
microscopy. This code was created by the David Grier group at NYU.

The best way to understand how it works is to see the code in
action. The example folder contains a small video, and the code for
analyzing the video. Simply go to dhm/example/programs and run

IDL>details_hd5f

and

IDL>rundhmstream.pro


Dependencies
-----------------------------------
dhm needs the following other repositories on my github:
features
math
utility
video.

The video repository needs mplayer and mencoder installed.  IDL also
has it's on version of certain codecs. These need to be removed in
order for the video code to work. The following files in (on my Ubuntu
setup) /usr/local/exelis/idl83/bin/bin.linux.x86_64 should be deleted:
libavcodec.so.54
libavcodec.so
libavformat.so.54
libavformat.so
libavutil.52
libavutil.so
libswscale.so.2
libswscale.so