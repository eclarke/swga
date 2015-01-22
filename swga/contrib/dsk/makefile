SVNREV = -DSVN_REV=1.6906
CC=g++
PROGNAME=dsk
CFLAGS+= -O4
LDFLAGS+=  -lz $(SVNREV)
SRC=minia/Pool.cpp minia/Bank.cpp minia/Bloom.cpp minia/Hash16.cpp minia/Kmer.cpp minia/Set.cpp minia/Utils.cpp minia/SortingCount.cpp  minia/LinearCounter.cpp minia/OAHash.cpp minia/LargeInt.cpp
EXEC=$(PROGNAME)
OBJ= $(SRC:.cpp=.o)
all: $(EXEC) parse_results

ifeq ($(prof),1)
 CFLAGS+=-pg
endif
ifeq ($(deb),1)
 #CFLAGS+=-O0 -DASSERTS -g -fstack-check # stack overflow protection
 CFLAGS+= -DASSERTS -g
endif


ifeq ($(omp),1)
 CFLAGS=-O4 -fopenmp -DOMP=1
endif

ifeq ($(osx),1)
 CFLAGS+= -DOSX=1
endif


k := 0$(k) # dummy k if not specified 
K_BELOW_32 := $(shell echo $(k)\<=32 | bc)
K_BELOW_64 := $(shell echo $(k)\<=64 | bc)
ARCH := $(shell getconf LONG_BIT) # detects sizeof(int)
USING_UINT128 := 0
largeintlib := 0

ifeq ($(K_BELOW_32),0)

    # use uint128 when 32<k<=64 and 64-bit architecture
    ifeq ($(K_BELOW_64),1)
        ifeq ($(strip $(ARCH)),64)
            CFLAGS += -Dkmer_type=__uint128_t
            USING_UINT128 := 1
        endif
    endif
    
    # use a bigint library otherwise
    ifeq ($(USING_UINT128),0)
        largeintlib := largeint#ttmath
    endif
endif

# ttmath (now, largeint) is used when you type "make k=[kmer size]" with a kmer size longer than supported integer type,
ifeq ($(largeintlib),ttmath)
    KMER_PRECISION := $(shell echo \($(k)+15\)/16 | bc)
endif
ifeq ($(largeintlib),largeint)
    KMER_PRECISION := $(shell echo \($(k)+31\)/32 | bc)
endif

ifneq ($(largeintlib),0)
    CFLAGS += -D_$(largeintlib) -DKMER_PRECISION=$(KMER_PRECISION)
endif

$(PROGNAME): clean $(OBJ) main.cpp
	$(CC) -o $@ $(OBJ) main.cpp $(CFLAGS) $(LDFLAGS)

parse_results: clean parse_results.cpp
	$(CC) -o $@ parse_results.cpp minia/Bank.cpp minia/Kmer.cpp minia/LargeInt.cpp $(CFLAGS) $(LDFLAGS)

%.o: %.cpp %.h
	$(CC) -o $@ -c $< $(CFLAGS)


%.o: %.c %.h 
	$(CC) -o $@ -c $< $(CFLAGS)
    
clean:
	rm -rf *.o minia/*.o
