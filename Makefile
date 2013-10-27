all:
	g++ program3.cpp  -DGL_GLEXT_PROTOTYPES -lGL -lglut -lGLU -lm
clean:
	rm *.out

