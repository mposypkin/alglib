OPTS = -O3 -I ../src
LIBS = ../src/alglib.a

all: tsph.exe tbc.exe
clean: 
	rm -f *.o *.exe

#tsph.exe:
#	c++ -O3 -I ../src -o tsph.exe tsph.cpp ../src/alglib.a
	
.o.exe:
	c++ $(OPTS) -o $@ $< $(LIBS) 

.cpp.o:
	c++ $(OPTS) -c $<


.SUFFIXES:
.SUFFIXES: .o .a .cpp .exe
