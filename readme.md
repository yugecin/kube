This is the source of an executable animation, released at the TRSAC 2022
demoparty.

![GIF render](kube150x150@20fps.gif?raw=true "GIF render")

also read 'kube.nfo'

https://demozoo.org/productions/314185/

The resulting animation can be seen in the 'kube.webm' file, or you can download
the release and run the exe yourself (beware - performance may be really bad).

The result can also be seen on ShaderToy: https://www.shadertoy.com/view/NtVfzd,
some shader code was edited for this to work - see the 'shadertoy' branch.

There are also rendered gifs available in the github release

note: this was made with a deadline and little time to that deadline, no  
effort was spent on cleanliness. expect chaos.

```
Tool credits
------------

| Crinkler - https://github.com/runestubbe/Crinkler
|
| "Crinkler is a compressing linker for Windows, specifically targeted
| towards executables with a size of just a few kilobytes. Its main purpose
| is as a tool for producing small demoscene productions."

File breakdown
--------------
Crinkler/
	contains crinkler
wip/
	contains random screenshots during the process of making this
c.c
	source of the program, copied and adapted from 'b.c' of my metro demo,
	see https://github.com/yugecin/metro/, so acknowledgements go out
	to auld's "chocolux" demo, iq's 4k framework, and sleep deprivation
frag.glsl
	source of the shader
frag.glsl.c
	frag.glsl but exported to be a c string, removing debug things,
	generated by https://github.com/yugecin/ShaderThing2
kube-entry-screenshot.png
	the image used to submit the entry in the partysystem
kube.nfo
	nfo file, basically a readme file that goes with the release,
	go read it
kube.webm
	a recording rendered to webm format
kube150x150@20fps.gif
	a recording rendered to gif format, see the release on GitHub for
	more sizes with different fps values
makedebug
	shell "script" to make a quick (aka debug) build
	unlike with my 'metro' demo, there is no 'makerelease' this time,
	this script was make to make the final released executables


a note on file modified dates: my laptop didn't have time correctly set up,
so most files were modified later than what it may look like, it was closer
to the deadline of 17:30


How 2 build the executable
--------------------------

Asumming a non-cmd terminal (git-bash or similar),
run ". makedebug", though you may need to edit the file to change the LIBPATH.

requires gcc

How frag.glsl.c is built
------------------------

Run https://github.com/yugecin/ShaderThing2 and open frag.glsl and press the
"export" button (in the main window - the one with the shader).

(Note that the tool might complain about a missing frag.glsl.pos file, ignore
it)
```
