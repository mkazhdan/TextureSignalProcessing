TARGET= LineIntegralConvolution
SOURCE= LineIntegralConvolution.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fpermissive -fopenmp -Wno-deprecated -Wno-unused-result -Wno-format -msse2 -std=c++17
	LFLAGS += -lgomp -lz -lGL -lGLU -lglut -lGLEW -lpng -ljpeg
	CC=gcc
	CXX=g++
else
	CFLAGS += -fpermissive -Wno-deprecated -Wno-unused-result -Wno-format -msse2 -std=c++17
	LFLAGS += -lgomp -lz -lGL -lGLU -lglut -lGLEW -lpng -ljpeg
	CC=clang
	CXX=clang++
endif

CFLAGS_DEBUG = -DDEBUG -g3
LFLAGS_DEBUG =

CFLAGS_RELEASE = -O3 -DRELEASE -funroll-loops -DNDEBUG
LFLAGS_RELEASE = -O3 

SRC = ./
BIN = ./../Bin/Linux/
BIN_O = ./
INCLUDE = /usr/include/ -I./../include/ -I./../third_party/include/ -I.

CC=gcc
CXX=g++
MD=mkdir

OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SOURCE))))

all: CFLAGS += $(CFLAGS_RELEASE)
all: LFLAGS += $(LFLAGS_RELEASE)
all: $(BIN)
all: $(BIN)$(TARGET)

debug: CFLAGS += $(CFLAGS_DEBUG)
debug: LFLAGS += $(LFLAGS_DEBUG)
debug: $(BIN)
debug: $(BIN)$(TARGET)

clean:
	rm -f $(BIN)$(TARGET)
	rm -f $(OBJECTS)

$(BIN):
	$(MD) -p $(BIN)

$(BIN)$(TARGET): $(OBJECTS)
	$(CXX) -o $@ $(OBJECTS) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.c
	$(CC) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<


