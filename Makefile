CXX      = clang++
CXXFLAGS = -std=c++11 -Wall -Wextra -I./include -DACCELERATE_NEW_LAPACK
LDFLAGS  = -framework Accelerate

SRC_DIR    = src
OBJ_DIR    = obj
BIN_DIR    = bin

SRCS        = $(wildcard $(SRC_DIR)/*.cpp)
TEST_SRC    = $(SRC_DIR)/test.cpp
TEST_SVD_SRC= $(SRC_DIR)/test_svd.cpp
TEST_OBJ    = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(TEST_SRC))
TEST_SVD_OBJ= $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(TEST_SVD_SRC))
LIB_OBJS    = $(filter-out $(TEST_OBJ) $(TEST_SVD_OBJ),$(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS)))

TARGET         = $(BIN_DIR)/sparse_matrix
TEST_TARGET    = $(BIN_DIR)/test
TEST_SVD_TARGET= $(BIN_DIR)/test_svd

.PHONY: all clean directories test test_svd

all: directories $(TARGET)

test: directories $(TEST_TARGET)
	$(TEST_TARGET)

test_svd: directories $(TEST_SVD_TARGET)
	$(TEST_SVD_TARGET)

directories:
	@mkdir -p $(OBJ_DIR) $(BIN_DIR)

$(TARGET): $(LIB_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TEST_TARGET): $(TEST_OBJ) $(LIB_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(TEST_SVD_TARGET): $(TEST_SVD_OBJ) $(LIB_OBJS)
	$(CXX) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
