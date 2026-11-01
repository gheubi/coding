# Polar Coding and Distribution Matching in Julia

This repository contains a Julia-based implementation for simulating and analyzing polar coding schemes, dynamic constraints, channel models, and distribution matching techniques.

The codebase is organized around a set of reusable utilities and simulation pipelines for:

- polar transforms and coding structures
- AWGN / BEC / BSC channel models
- CRC-constrained and dynamic frozen constraints
- successive-cancellation decoding logic
- lossy source coding simulations
- distribution matching and shaping experiments

## Repository layout

```text
.
├── conf.jl                    # Shared constants, tags, and configuration
├── data/                      # Simulation datasets and generated artifacts
├── src/
│   ├── main.jl                # Main entry point and project includes
│   ├── channel/               # Channel models and simulations
│   ├── channel_coding/        # Channel coding implementations and decoders
│   ├── distribution_matching/ # Distribution matching logic and simulations
│   ├── source_coding/         # Lossy source coding implementations
│   └── utils/                 # Generic helpers, entropy, parser, DB, CRC, etc.
├── test/
│   ├── polar_transform_test.jl
│   └── test_crc_constraint.jl
└── README.md
```

## Requirements

- Julia 1.12+ is recommended
- No package manager project file is included in this repository; the scripts are run directly from the workspace root

## Quick start

From the repository root, load the main file:

```bash
julia -e 'include("src/main.jl")'
```

This entry point includes the project modules and initializes the simulation utilities used across the codebase.

## Run tests

The project includes a test for the polar transform implementation:

```bash
julia test/polar_transform_test.jl
```

This test verifies bit-reversal behavior and the correctness of the polar transform against a matrix construction.

## Notes

- The repository includes a large number of configurable simulation parameters in `src/main.jl` and `conf.jl`.
- Channel and coding experiments are organized by module, allowing for targeted simulations and performance evaluation.
- The `data/` directory contains polar-code-related result files that can be used for comparison or further analysis.

## Typical use

This project is intended for research and simulation work around:

- polar code design and decoding
- constrained coding and dynamic frozen bits
- reliability analysis under noisy channels
- distribution matching and shaping for communication systems

## Development notes

The code is not currently packaged as a formal Julia package. Instead, it is structured as a lightweight simulation workspace where files are included directly from the `src/` tree.
