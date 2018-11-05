docker run 	--rm -it 	^
		-v %CD%\work:/home/pytfa/work 	^
		-v %CD%\..\..\pytfa:/src/pytfa		^
		-v %CD%\..:/src/etfl		^
		etfl-docker %*
