.PHONY: all clean build test linear_advection

# Compiler and flags
CXX := g++
CXXFLAGS := -Wall -Wextra -O2 -std=c++20 -Iinc
LDFLAGS :=

# Directories
SRC_DIR := src
BUILD_DIR := build
BIN_DIR := bin

# Targets
TEST_TARGET := $(BIN_DIR)/module_tests
LIN_ADV_TARGET := $(BIN_DIR)/linear_advection

# Default target
all: $(TEST_TARGET) $(LIN_ADV_TARGET)

# Create directories
$(BUILD_DIR) $(BIN_DIR):
	@mkdir -p $@

# Compile source file: module_tests
$(BUILD_DIR)/module_tests.o: $(SRC_DIR)/module_tests.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile source file: linear_advection
$(BUILD_DIR)/linear_advection.o: $(SRC_DIR)/linear_advection.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link executable module_tests
$(TEST_TARGET): $(BUILD_DIR)/module_tests.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# Link executable linear_advection
$(LIN_ADV_TARGET): $(BUILD_DIR)/linear_advection.o | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@

# Run the program
test: $(TEST_TARGET)
	./$(TEST_TARGET)

linear_advection: $(LIN_ADV_TARGET)
	./$(LIN_ADV_TARGET)

# Clean build artifacts
clean:
	@rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Cleaned build artifacts"

# Help target
help:
	@echo "Available targets:"
	@echo "  all   - Build both test and linear_advection target"
	@echo "  test  - Build and run module_tests.cpp"
	@echo "	 linear_advection - Build and run linear_advection.cpp"
	@echo "  clean - Remove build artifacts"
	@echo "  help  - Show this help message"
