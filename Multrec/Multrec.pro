QT       -= core
QT       -= gui

TARGET = Multrec
CONFIG   += console c++11
CONFIG   -= app_bundle

QMAKE_CXXFLAGS += -std=c++0x

TEMPLATE = app

SOURCES += main.cpp \
    trees/genespeciestreeutil.cpp \
    trees/newicklex.cpp \
    trees/node.cpp \
    trees/treeinfo.cpp \
    trees/treeiterator.cpp \
    multigenereconciler.cpp

HEADERS += \
    trees/genespeciestreeutil.h \
    trees/newicklex.h \
    trees/node.h \
    trees/treeinfo.h \
    trees/treeiterator.h \
    div/define.h \
    div/tinydir.h \
    div/util.h \
    multigenereconciler.h
