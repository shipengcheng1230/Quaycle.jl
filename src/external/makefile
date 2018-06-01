# Build Okada's code for calculating deformation due to a fault model
#
CC = gfortran
CFLAGS = -fPIC -w -O2 -g
LDFLAGS = -shared

SRCS = dc3d.f okada.f90
OBJS = $(SRCS:.c=.o)

TARGET = dc3d.so

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $(TARGET) $(OBJS)
