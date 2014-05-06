all: clean fractal.so

fractal.so: fractal.os	
	g++ -o fractal.so -shared fractal.os -L$(ARNOLD_PATH)/bin -lai

fractal.os: fractal.cpp
	g++ -o fractal.os -c -fPIC -D_LINUX -I$(ARNOLD_PATH)/include fractal.cpp

clean:
	rm -f fractal.os
	rm -f fractal.so

