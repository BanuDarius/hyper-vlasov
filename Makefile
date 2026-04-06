CC = gcc
CFLAGS = -Iinclude -fopenmp -flto -O3 -march=native -MMD -MP -s
LDLIBS = -lm -lgsl -lgslcblas

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
OUTPUT_DIR = output

TARGET = $(BIN_DIR)/hyper_vlasov

SRCS = $(SRC_DIR)/hyper_vlasov.c $(SRC_DIR)/init.c $(SRC_DIR)/tools.c $(SRC_DIR)/physics.c  $(SRC_DIR)/fit_algorithm.c

OBJS = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRCS))

all: $(TARGET) output_dirs

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(BUILD_DIR):
	mkdir -p $@

output_dirs:
	mkdir -p $(OUTPUT_DIR)

clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR) $(OUTPUT_DIR)

.PHONY: all clean