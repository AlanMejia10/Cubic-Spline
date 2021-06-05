TEMPLATE = app
CONFIG += console c++11
LIBS += -O2 -larmadillo
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp


QMAKE_CXXFLAGS += -Wno-padded
