# include conf/platform.mk

PYTHON=python2.7

SRCS = cexp.cpp

CFLAGS += -std=c++11 -Wall -Wno-conversion -Wno-unused-variable -Wno-return-type -O3 -ffast-math -fPIC
# CFLAGS += -static-libstdc++ -static-libgcc
# CFLAGS += -std=c++11 -Wall -Wno-conversion -Wno-unused-variable -Wno-return-type -O0 -g -Wfatal-errors -fPIC
LIBS += -lm
INCLUDES= 

OBJS = $(SRCS:.cpp=.o)

python: libcexp.so
	$(PYTHON) setup.py build_ext -i

.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<  -o $@

libcexp.so: $(OBJS)
	$(CXX) $(CFLAGS) $(INCLUDES) -o libcexp.so $(OBJS) $(LFLAGS) $(LIBS) -shared

clean:
	$(RM) order/*.o *.o *.so *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

headers:
	cd include && ../parse_maple.py 2


