TEMPLATE = app
CONFIG += console c++17
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
        src/functions.cpp \
        src/membrane_load_distribution.cpp \
        src/string_load_distribution.cpp \
        src/tests.cpp

HEADERS += \
    include/functions.h \
    include/log_duration.h \
    include/membrane_load_distribution.h \
    include/string_load_distribution.h \
    include/test_framework.h \
    include/tests.h \
    include/write.h
