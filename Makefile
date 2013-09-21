include base.mk

HADOOP_INSTALL_CAND = $(wildcard /home/kzhai/Public/hadoop-0.20.2 /usr/local/hadoop /usr/src/hadoop-0.20 /fs/cliplab/software/hadoop-0.20.2)
ifeq ($(HADOOP_INSTALL_CAND),)
    HADOOP_INSTALL = /usr/lib/hadoop
else
    HADOOP_INSTALL = $(HADOOP_INSTALL_CAND)
endif

PLATFORM = Linux-amd64-64
# If we're no longer using unordered_map, then we can remove last flag
CPPFLAGS = `gsl-config --cflags` -I../../../
#  file was built for unsupported file format
PIPESFLAGS = $(CPPFLAGS) -I$(HADOOP_INSTALL)/c++/$(PLATFORM)/include
#PIPESFLAGS = -m32 -I$(HADOOP_INSTALL)/c++/$(PLATFORM)/include

HADOOPLIBS = -L$(HADOOP_INSTALL)/c++/$(PLATFORM)/lib -lhadooppipes \
		-lhadooputils -lssl -lpthread -lcrypto -ldl
LIBS = $(LIBDIRS) $(HADOOPLIBS)

VOCAB_SIZE=5000

PROTO_OBJ = syntop_parameters.pb.o

$(PROTO_OBJ): proto/syntop_parameters.proto
	protoc $< -I=proto/ --cpp_out=src/ --python_out=.
	$(GPP) $(CPPFLAGS) $(INCLUDEDIRS) -Wall $(LIBDIRS) -c src/syntop_parameters.pb.cc

deploy_data: ../../data/wackypedia/hadoop/$(VOCAB_SIZE).0.prs
	hadoop dfs -rm wackypedia/$(VOCAB_SIZE).0.prs
	hadoop dfs -put ../../data/wackypedia/hadoop/$(VOCAB_SIZE).0.prs wackypedia

deploy_hdfs: wordcount
	hadoop dfs -rmr dft1-out
	hadoop dfs -rm bin/wordcount
	hadoop dfs -put wordcount bin/wordcount

LOCAL_CPP = src/util.cpp src/lens_doc.cpp src/gradient.cpp src/node.cpp src/vectorops.cpp src/document_mapper.cpp src/syntop_standalone_reducer.cpp src/syntop_standalone_mapper.cpp
LOCAL_OBJ = $(notdir $(LOCAL_CPP:.cpp=.o))
LOCAL_H = $(LOCAL_CPP:.cpp=.h)
LOCAL_DEP=$(wildcard src/*.*)

$(LOCAL_OBJ): $(LOCAL_CPP)
	$(GPP) $(PIPESFLAGS) $(CFLAGS) $(INCLUDEDIRS) $(LIBDIRS) -c $(LOCAL_CPP)

../../data/wackypedia/hadoop/$(VOCAB_SIZE).0.prs: scripts/wacky_reducer.py parse_reader.py $(WACKYPEDIA_NUMERIC)
	rm -f ../../data/wackypedia/hadoop/$(VOCAB_SIZE).*
	mkdir -p /tmp/`whoami`/syntop/wackypedia
	$(PYTHON_COMMAND) scripts/wacky_reducer.py --vocab_size=$(VOCAB_SIZE) --vocab_source="../../data/wackypedia/numeric/wpdia*.index" --output_path=/tmp/`whoami`/syntop/wackypedia/$(VOCAB_SIZE)
	mkdir -p ../../data/wackypedia/hadoop
	mv /tmp/`whoami`/syntop/wackypedia/$(VOCAB_SIZE).* ../../data/wackypedia/hadoop/

syntop: src/syntop.cpp $(PROTO_OBJ) $(LOCAL_OBJ) $(UTIL_OBJ)
	cpplint.py $(LINT_OPTIONS) $(LOCAL_H) $(LOCAL_CPP) src/syntop_mapper.h src/syntop_mapper.cpp src/syntop_reducer.h src/syntop_reducer.cpp
	$(GPP) $(PIPESFLAGS) $(INCLUDEDIRS) $(LOCAL_OBJ) $(PROTO_OBJ) $(UTIL_OBJ) src/syntop_reducer.cpp src/syntop_mapper.cpp $< $(LIBDIRS) -L$(HADOOP_INSTALL)/c++/$(PLATFORM)/lib $(HADOOPLIBS) -g -O2 -static -o $@

simulation: src/simulation.cpp $(PROTO_OBJ) $(LOCAL_OBJ) $(UTIL_OBJ)
	cpplint.py $(LINT_OPTIONS) $(LOCAL_H) $(LOCAL_CPP)
	$(GPP) $(PIPESFLAGS) $(INCLUDEDIRS) $(LOCAL_OBJ) $(PROTO_OBJ) $(UTIL_OBJ) $< $(LIBDIRS) -g -DEBUG  --save-temps -ggdb -o $@

syntop_test: src/test.cpp $(PROTO_OBJ) $(LOCAL_OBJ) $(UTIL_OBJ)
	cpplint.py $(LINT_OPTIONS) $(LOCAL_H) $(LOCAL_CPP)
	$(GPP) $(PIPESFLAGS) $(INCLUDEDIRS) $(LOCAL_OBJ) $(PROTO_OBJ) $(UTIL_OBJ) $< $(LIBDIRS) -L$(HADOOP_INSTALL)/c++/$(PLATFORM)/lib $(HADOOPLIBS) -o $@

#wordcount: src/wordcount.cpp $(PROTO_OBJ) $(LOCAL_OBJ) $(UTIL_OBJ)
#	cpplint.py $(LINT_OPTIONS) src/*.cpp src/*.h
#	$(GPP) $(PIPESFLAGS) $(INCLUDEDIRS) $(LOCAL_OBJ) $(PROTO_OBJ) $(UTIL_OBJ) $< $(LIBDIRS) -L$(HADOOP_INSTALL)/c++/$(PLATFORM)/lib -lhadooppipes -lhadooputils -lpthread -g -O2 -o $@

#batch_syntop: src/batch_syntop.cpp $(PROTO_OBJ) $(LOCAL_OBJ) $(UTIL_OBJ)
#	$(GPP) $(PIPESFLAGS) $(INCLUDEDIRS) $(LOCAL_OBJ) $(PROTO_OBJ) $(UTIL_OBJ) $< $(LIBDIRS) -g -DEBUG --save-temps -ggdb -o $@

#batch_wordcount: src/batch_wordcount.cpp $(LOCAL_OBJ) $(UTIL_OBJ)
#	$(GPP) $(PIPESFLAGS) $(INCLUDEDIRS) $(LOCAL_OBJ) $(UTIL_OBJ) $< $(LIBDIRS) -g -DEBUG --save-temps -ggdb -o $@
