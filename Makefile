CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -I./include
LDFLAGS = 

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
TEST_SRC = $(SRC_DIR)/test.cpp
TEST_OBJ = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(TEST_SRC))
LIB_OBJS = $(filter-out $(TEST_OBJ),$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS)))

TARGET = $(BIN_DIR)/sparse_matrix
TEST_TARGET = $(BIN_DIR)/test

.PHONY: all clean directories test

all: directories $(TARGET)

test: directories $(TEST_TARGET)
	$(TEST_TARGET)

directories:
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(BIN_DIR)

$(TARGET): $(LIB_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(TEST_TARGET): $(TEST_OBJ) $(LIB_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) 