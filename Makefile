CC = g++
OPT_FLAG = -O3
CFLAGS = $(OPT_FLAG) -march=native -Iinclude -fopenmp -flto -MMD -MP -g -Wall -Wextra -Wshadow
LDLIBS = -lm -lgsl

SRC_DIR = src
BUILD_DIR = build
BIN_DIR = bin
INPUT_DIR = input
OUTPUT_DIR = output
OUTPUT_IMAGE = output-image

TARGET = $(BIN_DIR)/hyper_vlasov

SRCS = $(SRC_DIR)/hyper_vlasov.cpp
OBJS = $(BUILD_DIR)/hyper_vlasov.o

all: output-dirs $(TARGET)

fast: OPT_FLAG = -Ofast
fast: all

$(TARGET): $(OBJS) | $(BIN_DIR)
	@$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)
	$(info Linked $@.)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@$(CC) $(CFLAGS) -c $< -o $@
	$(info Compiled $@.)

$(BIN_DIR) $(BUILD_DIR):
	@mkdir -p $@
	$(info Created $@ directory.)

output-dirs:
	@mkdir -p $(OUTPUT_DIR) $(OUTPUT_IMAGE) $(INPUT_DIR) $(BUILD_DIR $(BIN_DIR)):
	$(info Created output directories.)

clean:
	@rm -rf $(BUILD_DIR) $(BIN_DIR) $(OUTPUT_DIR) $(OUTPUT_IMAGE) $(INPUT_DIR)
	$(info Removed output files.)

-include $(OBJS:.o=.d)

.PHONY: all clean output_dirs fast