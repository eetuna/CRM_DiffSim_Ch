# Makefile for convenient commands

.PHONY: help setup data train evaluate test pipeline clean build-cpp gpu check

help:
	@echo "MRI Catheter ML/RL Makefile"
	@echo ""
	@echo "Main commands:"
	@echo "  make setup      - Complete setup (build bindings, install package)"
	@echo "  make data       - Generate training data"
	@echo "  make train      - Train all models"
	@echo "  make evaluate   - Evaluate trained agent"
	@echo "  make test       - Run tests"
	@echo "  make pipeline   - Run complete pipeline (data + train + eval + test)"
	@echo "  make clean      - Remove all generated files"
	@echo ""
	@echo "Development commands:"
	@echo "  make build-cpp  - Build C++ bindings only"
	@echo "  make gpu        - Enable GPU support (install CUDA PyTorch)"
	@echo "  make check      - Check system requirements"
	@echo ""
	@echo "Environment variables:"
	@echo "  CRM_PATH       - Path to CRM code (default: ../CRM)"
	@echo "  CONFIG_DIM     - Configuration dimension (default: 6)"
	@echo "  STATE_DIM      - State dimension (default: 24)"
	@echo ""
	@echo "Examples:"
	@echo "  CRM_PATH=/path/to/CRM make setup"
	@echo "  make build-cpp --clean      # Clean rebuild of bindings"
	@echo "  make gpu                    # Enable GPU acceleration"
	@echo "  make check                  # Verify system setup"
	@echo "  make pipeline               # Run everything"

setup:
	@chmod +x scripts/*.sh
	@./scripts/setup.sh

data:
	@./scripts/generate_data.sh

train:
	@./scripts/train_all.sh

evaluate:
	@./scripts/evaluate.sh

test:
	@./scripts/test.sh

pipeline:
	@./scripts/run_pipeline.sh

clean:
	@./scripts/clean.sh

build-cpp:
	@./scripts/build_cpp.sh $(filter-out $@,$(MAKECMDGOALS))

gpu:
	@./scripts/enable_gpu.sh

check:
	@./scripts/check_system.sh

# Allow passing --clean to build-cpp
%:
	@:
