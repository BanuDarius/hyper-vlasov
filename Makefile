CC = g++
OPT_FLAG = -O3
CFLAGS = -std=c++20 -Iinclude -fopenmp -flto $(OPT_FLAG) -march=native -MMD -MP -g -Wall -Wextra -Wshadow
LDLIBS = -lm -lgsl -lgslcblas

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
OUTPUT_DIR = output

TARGET = $(BIN_DIR)/hyper_vlasov

SRCS = $(SRC_DIR)/hyper_vlasov.cpp
OBJS = $(BUILD_DIR)/hyper_vlasov.o

all: output_dirs $(TARGET)

fast: OPT_FLAG = -Ofast
fast: all

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

output_dirs:
	mkdir -p $(OUTPUT_DIR)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) $(OUTPUT_DIR)

-include $(OBJS:.o=.d)

.PHONY: all clean output_dirs fast